import json
import os
import subprocess
import sys

import setuptools
import skbuild
from setuptools_scm import get_version

version = get_version(relative_to=__file__, write_to="python/resdata/version.py")


# Corporate networks tend to be behind a proxy server with their own non-public
# SSL certificates. Conan keeps its own certificates, whose path we can override
if "CONAN_CACERT_PATH" not in os.environ:
    # Look for a RHEL-compatible system-wide file
    for file_ in ("/etc/pki/tls/cert.pem",):
        if not os.path.isfile(file_):
            continue
        os.environ["CONAN_CACERT_PATH"] = file_
        break


def get_skbuild_dir():
    """Get the scikit-build cmake build directory path."""
    import platform as plat

    python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    return os.path.abspath(
        os.path.join(
            "_skbuild",
            f"{plat.system().lower()}-{plat.machine()}-{python_version}",
            "cmake-build",
        )
    )


def run_conan_install():
    """Run conan install to generate CMake presets and toolchain file for Conan 2."""
    import shutil

    conan_exe = shutil.which("conan")
    subprocess.run(
        [conan_exe, "profile", "detect", "--force"],
        check=False,  # Ignore if profile already exists
    )

    skbuild_dir = get_skbuild_dir()
    os.makedirs(skbuild_dir, exist_ok=True)

    subprocess.run(
        [
            conan_exe,
            "install",
            ".",
            f"--output-folder={skbuild_dir}",
            "--build=missing",
        ],
        check=True,
    )

    return skbuild_dir


def get_cmake_args_from_preset(skbuild_dir):
    """Extract CMake arguments from Conan-generated CMakePresets.json."""
    presets_file = os.path.join(skbuild_dir, "CMakePresets.json")
    cmake_args = []

    if os.path.exists(presets_file):
        with open(presets_file) as f:
            presets = json.load(f)

        if presets.get("configurePresets"):
            preset = presets["configurePresets"][0]

            if "toolchainFile" in preset:
                toolchain = preset["toolchainFile"]
                if not os.path.isabs(toolchain):
                    toolchain = os.path.join(skbuild_dir, toolchain)
                cmake_args.append(f"-DCMAKE_TOOLCHAIN_FILE='{toolchain}'")

            for key, value in preset.get("cacheVariables", {}).items():
                if isinstance(value, dict):
                    value = value.get("value", "")
                cmake_args.append(f"-D{key}={value}")

            for key, value in preset.get("environment", {}).items():
                if "$penv{" in value:
                    import re

                    value = re.sub(
                        r"\$penv\{(\w+)\}",
                        lambda m: os.environ.get(m.group(1), ""),
                        value,
                    )
                os.environ[key] = value

    return cmake_args


skbuild_dir = run_conan_install()
CMAKE_ARGS_FROM_PRESET = get_cmake_args_from_preset(skbuild_dir)


with open("README.md") as f:
    long_description = f.read()


def utility_wrappers():
    """
    Wrappers around resdata's "application" utilities. These are only supported on
    Linux at this time so only create the wrapper when on Linux.
    """
    if sys.platform != "linux":
        return []

    return [
        name + " = resdata.bin:main"
        for name in (
            "convert.x",  # deprecated
            "rd_pack.x",
            "rd_unpack.x",
            "grdecl_test.x",  # deprecated
            "kw_list.x",  # deprecated
            "load_test.x",  # deprecated
        )
    ] + ["summary.x = view_summary.__main__:main"]


skbuild.setup(
    name="resdata",
    author="Equinor ASA",
    author_email="fg_sib-scout@equinor.com",
    description="Package for reading and writing the fortran result files from reservoir simulators",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/equinor/resdata",
    packages=setuptools.find_packages(
        where="python",
        exclude=["*.tests", "*.tests.*", "tests.*", "tests", "ert.*", "ert"],
    ),
    package_dir={"": "python"},
    license="GPL-3.0",
    platforms="any",
    install_requires=[
        "cwrap",
        "numpy",
        "pandas",
        "natsort",
        "typing_extensions",
        "resfo-utilities>=0.4.0",
    ],
    setup_requires=["conan>=2"],
    entry_points={"console_scripts": utility_wrappers()},
    cmake_args=CMAKE_ARGS_FROM_PRESET
    + [
        "-DRD_VERSION=" + version,
        "-DBUILD_APPLICATIONS=" + ("ON" if sys.platform == "linux" else "OFF"),
        "-DBUILD_TESTS=OFF",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_INSTALL_BINDIR=python/resdata/.bin",
        "-DCMAKE_INSTALL_LIBDIR=python/resdata/.libs",
        "-DCMAKE_INSTALL_INCLUDEDIR=python/resdata/.include",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Programming Language :: Python :: 3.15",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries",
        "Topic :: Utilities",
    ],
    version=version,
)
