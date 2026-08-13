#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <htslib/sam.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct Config {
    int min_mapq;
    int min_baseq;
    int quality_radius;
    double quality_fraction;
    int min_clip_length;
    int min_indel_length;
    int min_aligned_length;
    int max_sa_nm;
    int discordant_min_distance;
    int max_insert_size;
    int merge_distance;
    bool allow_duplicates;
    bool include_supplementary;
    std::string library_orientation;
};

struct Clip {
    std::string chrom, side, sequence, read_name, strand, cigar, type;
    int64_t pos, length, reference_start, reference_end;
    int mapq;
    bool is_reverse;
};

struct Indel {
    std::string chrom, operation, sequence, read_name, cigar;
    int64_t start, end, length;
    int mapq;
};

struct Pair {
    std::string read_name, chrom, mate_chrom, orientation, cigar, reason;
    int64_t pos, mate_pos;
    int mapq;
    bool is_reverse, mate_is_reverse;
};

struct Split {
    std::string read_name, chrom, side, remote_chrom, remote_strand;
    std::string remote_cigar, orientation, cigar, sa_tag;
    int64_t pos, remote_pos;
    int remote_mapq, remote_nm, mapq;
};

struct Evidence {
    std::vector<Clip> clips;
    std::vector<Indel> indels;
    std::vector<Pair> pairs;
    std::vector<Split> splits;
};

struct BamReader {
    samFile* file = nullptr;
    sam_hdr_t* header = nullptr;
    hts_idx_t* index = nullptr;
    bam1_t* record = nullptr;

    explicit BamReader(const char* path) {
        file = sam_open(path, "rb");
        if (file) header = sam_hdr_read(file);
        if (header) index = sam_index_load(file, path);
        if (index) record = bam_init1();
    }

    ~BamReader() {
        if (record) bam_destroy1(record);
        if (index) hts_idx_destroy(index);
        if (header) sam_hdr_destroy(header);
        if (file) sam_close(file);
    }

    bool valid() const { return file && header && index && record; }
};

bool dict_long(PyObject* dict, const char* key, int& out) {
    PyObject* value = PyDict_GetItemString(dict, key);
    if (!value) {
        PyErr_Format(PyExc_KeyError, "Missing native evidence setting: %s", key);
        return false;
    }
    long parsed = PyLong_AsLong(value);
    if (parsed == -1 && PyErr_Occurred()) return false;
    out = static_cast<int>(parsed);
    return true;
}

bool dict_bool(PyObject* dict, const char* key, bool& out) {
    PyObject* value = PyDict_GetItemString(dict, key);
    if (!value) {
        PyErr_Format(PyExc_KeyError, "Missing native evidence setting: %s", key);
        return false;
    }
    int parsed = PyObject_IsTrue(value);
    if (parsed < 0) return false;
    out = parsed != 0;
    return true;
}

bool parse_config(PyObject* dict, int min_mapq, int min_baseq, Config& out) {
    if (!PyDict_Check(dict)) {
        PyErr_SetString(PyExc_TypeError, "config must be a dict");
        return false;
    }
    out.min_mapq = min_mapq;
    out.min_baseq = min_baseq;
    if (!dict_long(dict, "breakpoint_quality_window", out.quality_radius) ||
        !dict_long(dict, "min_clip_length", out.min_clip_length) ||
        !dict_long(dict, "min_indel_length", out.min_indel_length) ||
        !dict_long(dict, "min_aligned_length", out.min_aligned_length) ||
        !dict_long(dict, "max_sa_nm", out.max_sa_nm) ||
        !dict_long(dict, "discordant_min_distance", out.discordant_min_distance) ||
        !dict_long(dict, "max_insert_size", out.max_insert_size) ||
        !dict_long(dict, "merge_distance", out.merge_distance) ||
        !dict_bool(dict, "allow_duplicates", out.allow_duplicates) ||
        !dict_bool(dict, "include_supplementary", out.include_supplementary)) {
        return false;
    }
    PyObject* fraction = PyDict_GetItemString(dict, "min_breakpoint_quality_fraction");
    PyObject* orientation = PyDict_GetItemString(dict, "library_orientation");
    if (!fraction || !orientation) {
        PyErr_SetString(PyExc_KeyError, "Missing native evidence quality/orientation setting");
        return false;
    }
    out.quality_fraction = PyFloat_AsDouble(fraction);
    const char* text = PyUnicode_AsUTF8(orientation);
    if (PyErr_Occurred() || !text) return false;
    out.library_orientation = text;
    for (char& ch : out.library_orientation) {
        if (ch >= 'A' && ch <= 'Z') ch = static_cast<char>(ch - 'A' + 'a');
    }
    return true;
}

bool consumes_reference(int op) { return (bam_cigar_type(op) & 2) != 0; }
bool consumes_query(int op) { return (bam_cigar_type(op) & 1) != 0; }
bool aligned_query_op(int op) { return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF; }

int64_t aligned_query_length(const bam1_t* record) {
    const uint32_t* cigar = bam_get_cigar(record);
    int64_t length = 0;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        if (aligned_query_op(bam_cigar_op(cigar[i]))) length += bam_cigar_oplen(cigar[i]);
    }
    return length;
}

bool passes_read_filters(const bam1_t* record, const Config& config) {
    const uint16_t flag = record->core.flag;
    if ((flag & BAM_FUNMAP) || (flag & BAM_FSECONDARY)) return false;
    if ((flag & BAM_FSUPPLEMENTARY) && !config.include_supplementary) return false;
    if ((flag & BAM_FDUP) && !config.allow_duplicates) return false;
    if (record->core.n_cigar == 0 || record->core.qual < config.min_mapq) return false;
    return aligned_query_length(record) >= config.min_aligned_length;
}

std::string cigar_string(const bam1_t* record) {
    const uint32_t* cigar = bam_get_cigar(record);
    std::ostringstream out;
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        out << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return out.str().empty() ? "*" : out.str();
}

std::string query_sequence(const bam1_t* record) {
    const uint8_t* encoded = bam_get_seq(record);
    std::string sequence;
    sequence.resize(record->core.l_qseq);
    for (int i = 0; i < record->core.l_qseq; ++i) {
        sequence[i] = seq_nt16_str[bam_seqi(encoded, i)];
    }
    return sequence;
}

bool passes_breakpoint_quality(const bam1_t* record, int query_pos, const Config& config) {
    if (config.min_baseq <= 0) return true;
    const int radius = config.quality_radius;
    const int query_length = record->core.l_qseq;
    if (query_pos < radius || query_pos + radius > query_length || query_length == 0) return false;
    const uint8_t* qualities = bam_get_qual(record);
    if (!qualities || qualities[0] == 0xff) return false;
    int passing = 0;
    for (int i = query_pos - radius; i < query_pos + radius; ++i) {
        if (qualities[i] >= config.min_baseq) ++passing;
    }
    const int required = static_cast<int>(std::ceil(2.0 * radius * config.quality_fraction));
    return passing >= required;
}

struct ClipLengths {
    int left = 0, right = 0;
    char left_type = '\0', right_type = '\0';
};

ClipLengths clip_lengths(const std::vector<std::pair<int, int>>& cigar) {
    ClipLengths result;
    if (cigar.empty()) return result;
    const auto [first_op, first_length] = cigar.front();
    const auto [last_op, last_length] = cigar.back();
    if (first_op == BAM_CSOFT_CLIP || first_op == BAM_CHARD_CLIP) {
        result.left = first_length;
        result.left_type = first_op == BAM_CSOFT_CLIP ? 'S' : 'H';
    }
    if (last_op == BAM_CSOFT_CLIP || last_op == BAM_CHARD_CLIP) {
        result.right = last_length;
        result.right_type = last_op == BAM_CSOFT_CLIP ? 'S' : 'H';
    }
    return result;
}

std::vector<std::pair<int, int>> record_cigar(const bam1_t* record) {
    std::vector<std::pair<int, int>> result;
    result.reserve(record->core.n_cigar);
    const uint32_t* cigar = bam_get_cigar(record);
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        result.emplace_back(bam_cigar_op(cigar[i]), bam_cigar_oplen(cigar[i]));
    }
    return result;
}

int cigar_code(char op) {
    switch (op) {
        case 'M': return BAM_CMATCH;
        case 'I': return BAM_CINS;
        case 'D': return BAM_CDEL;
        case 'N': return BAM_CREF_SKIP;
        case 'S': return BAM_CSOFT_CLIP;
        case 'H': return BAM_CHARD_CLIP;
        case 'P': return BAM_CPAD;
        case '=': return BAM_CEQUAL;
        case 'X': return BAM_CDIFF;
        default: return -1;
    }
}

bool parse_cigar(const std::string& text, std::vector<std::pair<int, int>>& result) {
    int length = 0;
    bool have_digits = false;
    for (char ch : text) {
        if (ch >= '0' && ch <= '9') {
            length = length * 10 + (ch - '0');
            have_digits = true;
            continue;
        }
        const int op = cigar_code(ch);
        if (!have_digits || op < 0) return false;
        result.emplace_back(op, length);
        length = 0;
        have_digits = false;
    }
    return !have_digits && !result.empty();
}

struct Geometry {
    int64_t query_start, query_end, reference_start, reference_end;
};

Geometry alignment_geometry(int64_t reference_start,
                            const std::vector<std::pair<int, int>>& cigar,
                            char strand) {
    const ClipLengths clips = clip_lengths(cigar);
    int64_t query_span = 0;
    int64_t reference_span = 0;
    for (const auto [op, length] : cigar) {
        if (consumes_query(op) || op == BAM_CHARD_CLIP) query_span += length;
        if (consumes_reference(op)) reference_span += length;
    }
    const int64_t query_start = strand == '+' ? clips.left : clips.right;
    const int64_t query_end = strand == '+' ? query_span - clips.right : query_span - clips.left;
    return {query_start, query_end, reference_start, reference_start + reference_span};
}

std::pair<int64_t, std::string> breakpoint_for_query_edge(const Geometry& geometry,
                                                           char strand,
                                                           bool start_edge) {
    const bool uses_reference_start = start_edge == (strand == '+');
    return uses_reference_start
               ? std::make_pair(geometry.reference_start, std::string("left_clip"))
               : std::make_pair(geometry.reference_end, std::string("right_clip"));
}

void extract_record(const bam1_t* record,
                    const sam_hdr_t* header,
                    int64_t region_start,
                    int64_t region_end,
                    const Config& config,
                    Evidence& evidence,
                    std::unordered_set<std::string>& seen_pairs,
                    std::unordered_set<std::string>& seen_splits) {
    const std::string chrom = sam_hdr_tid2name(header, record->core.tid);
    const std::string read_name = bam_get_qname(record);
    const std::string cigar = cigar_string(record);
    const std::string sequence = query_sequence(record);
    const int64_t reference_start = record->core.pos;
    const int64_t reference_end = bam_endpos(record);
    const bool is_reverse = (record->core.flag & BAM_FREVERSE) != 0;
    const bool mate_is_reverse = (record->core.flag & BAM_FMREVERSE) != 0;
    const std::string strand = is_reverse ? "-" : "+";
    const auto cigar_ops = record_cigar(record);
    const ClipLengths clips = clip_lengths(cigar_ops);

    if (clips.left >= config.min_clip_length &&
        passes_breakpoint_quality(record,
                                  clips.left_type == 'S' ? clips.left : -1,
                                  config)) {
        evidence.clips.push_back({chrom,
                                  "left_clip",
                                  clips.left_type == 'S' && !sequence.empty()
                                      ? sequence.substr(0, clips.left)
                                      : "NA",
                                  read_name,
                                  strand,
                                  cigar,
                                  std::string(1, clips.left_type),
                                  reference_start,
                                  clips.left,
                                  reference_start,
                                  reference_end,
                                  record->core.qual,
                                  is_reverse});
    }
    const int right_boundary = clips.right_type == 'S'
                                   ? static_cast<int>(sequence.size()) - clips.right
                                   : -1;
    if (clips.right >= config.min_clip_length &&
        passes_breakpoint_quality(record, right_boundary, config)) {
        evidence.clips.push_back({chrom,
                                  "right_clip",
                                  clips.right_type == 'S' && !sequence.empty()
                                      ? sequence.substr(sequence.size() - clips.right)
                                      : "NA",
                                  read_name,
                                  strand,
                                  cigar,
                                  std::string(1, clips.right_type),
                                  reference_end,
                                  clips.right,
                                  reference_start,
                                  reference_end,
                                  record->core.qual,
                                  is_reverse});
    }

    int64_t ref_pos = reference_start;
    int query_pos = 0;
    for (const auto [op, length] : cigar_ops) {
        if (op == BAM_CINS && length >= config.min_indel_length &&
            passes_breakpoint_quality(record, query_pos, config) &&
            passes_breakpoint_quality(record, query_pos + length, config)) {
            const std::string inserted = sequence.empty()
                                             ? "NA"
                                             : sequence.substr(query_pos, length);
            evidence.indels.push_back({chrom,
                                       "INS",
                                       inserted.empty() ? "NA" : inserted,
                                       read_name,
                                       cigar,
                                       ref_pos,
                                       ref_pos,
                                       length,
                                       record->core.qual});
        } else if (op == BAM_CDEL && length >= config.min_indel_length &&
                   passes_breakpoint_quality(record, query_pos, config)) {
            evidence.indels.push_back({chrom,
                                       "DEL",
                                       "NA",
                                       read_name,
                                       cigar,
                                       ref_pos,
                                       ref_pos + length,
                                       length,
                                       record->core.qual});
        }
        if (consumes_reference(op)) ref_pos += length;
        if (consumes_query(op)) query_pos += length;
    }

    const uint16_t flag = record->core.flag;
    if ((flag & BAM_FPAIRED) && !(flag & BAM_FMUNMAP)) {
        const char* mate_name = record->core.mtid >= 0
                                    ? sam_hdr_tid2name(header, record->core.mtid)
                                    : nullptr;
        const std::string mate_chrom = mate_name ? mate_name : chrom;
        const int64_t mate_pos = record->core.mpos;
        std::vector<std::string> reasons;
        if (chrom != mate_chrom) {
            reasons.emplace_back("different_chrom");
        } else {
            if (std::llabs(mate_pos - reference_start) > config.discordant_min_distance)
                reasons.emplace_back("distant_mate");
            if (std::llabs(record->core.isize) > config.max_insert_size)
                reasons.emplace_back("large_insert");
        }
        const bool same_strand = is_reverse == mate_is_reverse;
        if ((config.library_orientation == "fr" || config.library_orientation == "rf") &&
            same_strand) {
            reasons.emplace_back("same_strand_pair");
        } else if ((config.library_orientation == "ff" ||
                    config.library_orientation == "rr") &&
                   !same_strand) {
            reasons.emplace_back("opposite_strand_pair");
        }
        const bool mate_outside = chrom == mate_chrom &&
                                  (mate_pos < region_start - config.merge_distance ||
                                   mate_pos > region_end + config.merge_distance);
        if (mate_outside && !reasons.empty())
            reasons.emplace_back("mate_outside_candidate_region");
        if (!reasons.empty() && seen_pairs.insert(read_name).second) {
            std::ostringstream reason;
            for (size_t i = 0; i < reasons.size(); ++i) {
                if (i) reason << ',';
                reason << reasons[i];
            }
            evidence.pairs.push_back({read_name,
                                      chrom,
                                      mate_chrom,
                                      std::string(is_reverse ? "-" : "+") +
                                          (mate_is_reverse ? "-" : "+"),
                                      cigar,
                                      reason.str(),
                                      reference_start,
                                      mate_pos,
                                      record->core.qual,
                                      is_reverse,
                                      mate_is_reverse});
        }
    }

    uint8_t* sa_value = bam_aux_get(record, "SA");
    const char* sa_text = sa_value ? bam_aux2Z(sa_value) : nullptr;
    if (!sa_text) return;
    const std::string sa_tag = sa_text;
    const char orientation_local = is_reverse ? '-' : '+';
    const Geometry local_geometry = alignment_geometry(reference_start, cigar_ops, orientation_local);

    size_t item_start = 0;
    while (item_start < sa_tag.size()) {
        const size_t item_end = sa_tag.find(';', item_start);
        const std::string item = sa_tag.substr(
            item_start, item_end == std::string::npos ? std::string::npos : item_end - item_start);
        item_start = item_end == std::string::npos ? sa_tag.size() : item_end + 1;
        if (item.empty()) continue;
        std::vector<std::string> fields;
        size_t field_start = 0;
        while (field_start <= item.size()) {
            const size_t comma = item.find(',', field_start);
            fields.push_back(item.substr(
                field_start, comma == std::string::npos ? std::string::npos : comma - field_start));
            if (comma == std::string::npos) break;
            field_start = comma + 1;
        }
        if (fields.size() < 6) continue;
        try {
            const std::string& remote_chrom = fields[0];
            const int64_t remote_start = std::stoll(fields[1]) - 1;
            const char remote_strand = fields[2].empty() ? '\0' : fields[2][0];
            const std::string& remote_cigar = fields[3];
            const int remote_mapq = std::stoi(fields[4]);
            const int remote_nm = std::stoi(fields[5]);
            std::vector<std::pair<int, int>> remote_ops;
            if ((remote_strand != '+' && remote_strand != '-') ||
                !parse_cigar(remote_cigar, remote_ops) ||
                remote_mapq < config.min_mapq || remote_nm > config.max_sa_nm) {
                continue;
            }
            const Geometry remote_geometry =
                alignment_geometry(remote_start, remote_ops, remote_strand);
            const bool local_first = local_geometry.query_start + local_geometry.query_end <=
                                     remote_geometry.query_start + remote_geometry.query_end;
            const bool local_start_edge = !local_first;
            const bool remote_start_edge = local_first;
            const auto [local_pos, side] = breakpoint_for_query_edge(
                local_geometry, orientation_local, local_start_edge);
            int query_boundary = -1;
            if (side == "left_clip" && clips.left_type == 'S')
                query_boundary = clips.left;
            else if (side == "right_clip" && clips.right_type == 'S')
                query_boundary = static_cast<int>(sequence.size()) - clips.right;
            if (!passes_breakpoint_quality(record, query_boundary, config)) continue;
            const auto [remote_pos, remote_side] = breakpoint_for_query_edge(
                remote_geometry, remote_strand, remote_start_edge);
            (void)remote_side;
            const std::string orientation = std::string(1, orientation_local) + remote_strand;
            const std::string key = read_name + '\x1f' + chrom + '\x1f' +
                                    std::to_string(local_pos) + '\x1f' + remote_chrom + '\x1f' +
                                    std::to_string(remote_pos) + '\x1f' + orientation;
            if (!seen_splits.insert(key).second) continue;
            evidence.splits.push_back({read_name,
                                       chrom,
                                       side,
                                       remote_chrom,
                                       std::string(1, remote_strand),
                                       remote_cigar,
                                       orientation,
                                       cigar,
                                       sa_tag,
                                       local_pos,
                                       remote_pos,
                                       remote_mapq,
                                       remote_nm,
                                       record->core.qual});
        } catch (const std::exception&) {
            continue;
        }
    }
}

PyObject* py_string(const std::string& value) {
    return PyUnicode_DecodeUTF8(value.data(), static_cast<Py_ssize_t>(value.size()), "strict");
}

PyObject* clip_tuple(const Clip& item) {
    PyObject* row = PyTuple_New(13);
    if (!row) return nullptr;
    PyTuple_SET_ITEM(row, 0, py_string(item.chrom));
    PyTuple_SET_ITEM(row, 1, PyLong_FromLongLong(item.pos));
    PyTuple_SET_ITEM(row, 2, py_string(item.side));
    PyTuple_SET_ITEM(row, 3, PyLong_FromLongLong(item.length));
    PyTuple_SET_ITEM(row, 4, py_string(item.sequence));
    PyTuple_SET_ITEM(row, 5, py_string(item.read_name));
    PyTuple_SET_ITEM(row, 6, py_string(item.strand));
    PyTuple_SET_ITEM(row, 7, PyLong_FromLong(item.mapq));
    PyTuple_SET_ITEM(row, 8, py_string(item.cigar));
    PyTuple_SET_ITEM(row, 9, PyBool_FromLong(item.is_reverse));
    PyTuple_SET_ITEM(row, 10, PyLong_FromLongLong(item.reference_start));
    PyTuple_SET_ITEM(row, 11, PyLong_FromLongLong(item.reference_end));
    PyTuple_SET_ITEM(row, 12, py_string(item.type));
    return row;
}

PyObject* indel_tuple(const Indel& item) {
    PyObject* row = PyTuple_New(9);
    if (!row) return nullptr;
    PyTuple_SET_ITEM(row, 0, py_string(item.chrom));
    PyTuple_SET_ITEM(row, 1, PyLong_FromLongLong(item.start));
    PyTuple_SET_ITEM(row, 2, PyLong_FromLongLong(item.end));
    PyTuple_SET_ITEM(row, 3, py_string(item.operation));
    PyTuple_SET_ITEM(row, 4, PyLong_FromLongLong(item.length));
    PyTuple_SET_ITEM(row, 5, py_string(item.sequence));
    PyTuple_SET_ITEM(row, 6, py_string(item.read_name));
    PyTuple_SET_ITEM(row, 7, PyLong_FromLong(item.mapq));
    PyTuple_SET_ITEM(row, 8, py_string(item.cigar));
    return row;
}

PyObject* pair_tuple(const Pair& item) {
    PyObject* row = PyTuple_New(11);
    if (!row) return nullptr;
    PyTuple_SET_ITEM(row, 0, py_string(item.read_name));
    PyTuple_SET_ITEM(row, 1, py_string(item.chrom));
    PyTuple_SET_ITEM(row, 2, PyLong_FromLongLong(item.pos));
    PyTuple_SET_ITEM(row, 3, py_string(item.mate_chrom));
    PyTuple_SET_ITEM(row, 4, PyLong_FromLongLong(item.mate_pos));
    PyTuple_SET_ITEM(row, 5, py_string(item.orientation));
    PyTuple_SET_ITEM(row, 6, PyLong_FromLong(item.mapq));
    PyTuple_SET_ITEM(row, 7, PyBool_FromLong(item.is_reverse));
    PyTuple_SET_ITEM(row, 8, PyBool_FromLong(item.mate_is_reverse));
    PyTuple_SET_ITEM(row, 9, py_string(item.cigar));
    PyTuple_SET_ITEM(row, 10, py_string(item.reason));
    return row;
}

PyObject* split_tuple(const Split& item) {
    PyObject* row = PyTuple_New(14);
    if (!row) return nullptr;
    PyTuple_SET_ITEM(row, 0, py_string(item.read_name));
    PyTuple_SET_ITEM(row, 1, py_string(item.chrom));
    PyTuple_SET_ITEM(row, 2, PyLong_FromLongLong(item.pos));
    PyTuple_SET_ITEM(row, 3, py_string(item.side));
    PyTuple_SET_ITEM(row, 4, py_string(item.remote_chrom));
    PyTuple_SET_ITEM(row, 5, PyLong_FromLongLong(item.remote_pos));
    PyTuple_SET_ITEM(row, 6, py_string(item.remote_strand));
    PyTuple_SET_ITEM(row, 7, py_string(item.remote_cigar));
    PyTuple_SET_ITEM(row, 8, PyLong_FromLong(item.remote_mapq));
    PyTuple_SET_ITEM(row, 9, PyLong_FromLong(item.remote_nm));
    PyTuple_SET_ITEM(row, 10, py_string(item.orientation));
    PyTuple_SET_ITEM(row, 11, PyLong_FromLong(item.mapq));
    PyTuple_SET_ITEM(row, 12, py_string(item.cigar));
    PyTuple_SET_ITEM(row, 13, py_string(item.sa_tag));
    return row;
}

template <typename T, PyObject* (*Convert)(const T&)>
PyObject* vector_list(const std::vector<T>& values) {
    PyObject* list = PyList_New(static_cast<Py_ssize_t>(values.size()));
    if (!list) return nullptr;
    for (Py_ssize_t i = 0; i < static_cast<Py_ssize_t>(values.size()); ++i) {
        PyObject* row = Convert(values[static_cast<size_t>(i)]);
        if (!row) {
            Py_DECREF(list);
            return nullptr;
        }
        PyList_SET_ITEM(list, i, row);
    }
    return list;
}

PyObject* evidence_dict(const Evidence& evidence) {
    PyObject* result = PyDict_New();
    PyObject* clips = vector_list<Clip, clip_tuple>(evidence.clips);
    PyObject* indels = vector_list<Indel, indel_tuple>(evidence.indels);
    PyObject* pairs = vector_list<Pair, pair_tuple>(evidence.pairs);
    PyObject* splits = vector_list<Split, split_tuple>(evidence.splits);
    if (!result || !clips || !indels || !pairs || !splits) {
        Py_XDECREF(result);
        Py_XDECREF(clips);
        Py_XDECREF(indels);
        Py_XDECREF(pairs);
        Py_XDECREF(splits);
        return nullptr;
    }
    PyDict_SetItemString(result, "clip_sites", clips);
    PyDict_SetItemString(result, "indels", indels);
    PyDict_SetItemString(result, "discordant_pairs", pairs);
    PyDict_SetItemString(result, "split_reads", splits);
    Py_DECREF(clips);
    Py_DECREF(indels);
    Py_DECREF(pairs);
    Py_DECREF(splits);
    return result;
}

PyObject* collect_regions(PyObject*, PyObject* args, PyObject* kwargs) {
    const char* bam_path = nullptr;
    PyObject* regions = nullptr;
    PyObject* config_dict = nullptr;
    int padding = 0;
    int min_mapq = 0;
    int min_baseq = 0;
    static const char* keywords[] = {
        "bam_path", "regions", "config", "padding", "min_mapq", "min_baseq", nullptr};
    if (!PyArg_ParseTupleAndKeywords(args,
                                     kwargs,
                                     "sOOiii:collect_regions",
                                     const_cast<char**>(keywords),
                                     &bam_path,
                                     &regions,
                                     &config_dict,
                                     &padding,
                                     &min_mapq,
                                     &min_baseq)) {
        return nullptr;
    }
    Config config;
    if (!parse_config(config_dict, min_mapq, min_baseq, config)) return nullptr;
    PyObject* region_list = PySequence_Fast(regions, "regions must be a sequence");
    if (!region_list) return nullptr;
    BamReader reader(bam_path);
    if (!reader.valid()) {
        Py_DECREF(region_list);
        PyErr_Format(PyExc_OSError, "Could not open indexed BAM with HTSlib: %s", bam_path);
        return nullptr;
    }
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(region_list);
    PyObject* result = PyList_New(count);
    if (!result) {
        Py_DECREF(region_list);
        return nullptr;
    }
    for (Py_ssize_t i = 0; i < count; ++i) {
        PyObject* region = PySequence_Fast(
            PySequence_Fast_GET_ITEM(region_list, i), "each region must be (chrom, start, end)");
        if (!region || PySequence_Fast_GET_SIZE(region) < 3) {
            Py_XDECREF(region);
            Py_DECREF(region_list);
            Py_DECREF(result);
            PyErr_SetString(PyExc_ValueError, "each region must be (chrom, start, end)");
            return nullptr;
        }
        const char* chrom = PyUnicode_AsUTF8(PySequence_Fast_GET_ITEM(region, 0));
        const int64_t region_start = PyLong_AsLongLong(PySequence_Fast_GET_ITEM(region, 1));
        const int64_t region_end = PyLong_AsLongLong(PySequence_Fast_GET_ITEM(region, 2));
        if (!chrom || PyErr_Occurred()) {
            Py_DECREF(region);
            Py_DECREF(region_list);
            Py_DECREF(result);
            return nullptr;
        }
        const int tid = sam_hdr_name2tid(reader.header, chrom);
        if (tid < 0) {
            Py_DECREF(region);
            Py_DECREF(region_list);
            Py_DECREF(result);
            PyErr_Format(PyExc_ValueError, "Invalid BAM contig: %s", chrom);
            return nullptr;
        }
        const int64_t fetch_start = std::max<int64_t>(0, region_start - padding);
        const int64_t fetch_end = region_end + padding;
        hts_itr_t* iterator = sam_itr_queryi(reader.index, tid, fetch_start, fetch_end);
        if (!iterator) {
            Py_DECREF(region);
            Py_DECREF(region_list);
            Py_DECREF(result);
            PyErr_Format(PyExc_RuntimeError, "Could not create BAM iterator for %s", chrom);
            return nullptr;
        }
        Evidence evidence;
        std::unordered_set<std::string> seen_pairs;
        std::unordered_set<std::string> seen_splits;
        int status = 0;
        while ((status = sam_itr_next(reader.file, iterator, reader.record)) >= 0) {
            if (passes_read_filters(reader.record, config)) {
                extract_record(reader.record,
                               reader.header,
                               region_start,
                               region_end,
                               config,
                               evidence,
                               seen_pairs,
                               seen_splits);
            }
        }
        hts_itr_destroy(iterator);
        Py_DECREF(region);
        if (status < -1) {
            Py_DECREF(region_list);
            Py_DECREF(result);
            PyErr_Format(PyExc_OSError, "HTSlib failed while reading %s:%lld-%lld",
                         chrom,
                         static_cast<long long>(fetch_start),
                         static_cast<long long>(fetch_end));
            return nullptr;
        }
        PyObject* item = evidence_dict(evidence);
        if (!item) {
            Py_DECREF(region_list);
            Py_DECREF(result);
            return nullptr;
        }
        PyList_SET_ITEM(result, i, item);
    }
    Py_DECREF(region_list);
    return result;
}

PyMethodDef methods[] = {
    {"collect_regions",
     reinterpret_cast<PyCFunction>(collect_regions),
     METH_VARARGS | METH_KEYWORDS,
     "Collect read evidence for a batch of candidate regions using HTSlib."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_evidence_native",
    "HTSlib-backed XXEJ evidence extraction.",
    -1,
    methods,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

}  // namespace

PyMODINIT_FUNC PyInit__evidence_native() { return PyModule_Create(&module); }
