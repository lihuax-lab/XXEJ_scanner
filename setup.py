from __future__ import annotations

import shlex
import subprocess

from setuptools import Extension, setup


def htslib_flags() -> tuple[list[str], list[str], list[str]]:
    try:
        flags = shlex.split(
            subprocess.check_output(
                ["pkg-config", "--cflags", "--libs", "htslib"], text=True
            )
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return [], [], ["hts"]
    include_dirs = [flag[2:] for flag in flags if flag.startswith("-I")]
    library_dirs = [flag[2:] for flag in flags if flag.startswith("-L")]
    libraries = [flag[2:] for flag in flags if flag.startswith("-l")]
    return include_dirs, library_dirs, libraries


include_dirs, library_dirs, libraries = htslib_flags()

setup(
    ext_modules=[
        Extension(
            "XXEJ_scanner._evidence_native",
            ["src/XXEJ_scanner/_evidence_native.cpp"],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            extra_compile_args=["-O3", "-std=c++17"],
            language="c++",
            optional=True,
        )
    ]
)
