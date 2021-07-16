import sys
import skbuild
import setuptools
from setuptools_scm import get_version


version = get_version(relative_to=__file__, write_to="python/ecl/version.py")


with open("README.md") as f:
    long_description = f.read()


def utility_wrappers():
    """
    Wrappers around ecl's "application" utilities. These are only supported on
    Linux at this time so only create the wrapper when on Linux.
    """
    if sys.platform != "linux":
        return []

    return [
        name + " = ecl.bin:main"
        for name in (
            "CF_dump",
            "convert.x",
            "ecl_pack.x",
            "ecl_unpack.x",
            "grdecl_grid",
            "grdecl_test.x",
            "grid_dump.x",
            "grid_dump_ascii.x",
            "grid_info.x",
            "kw_extract",
            "kw_list.x",
            "load_test.x",
            "make_grid",
            "ri_well_test",
            "segment_info",
            "select_test.x",
            "summary.x",
        )
    ]


skbuild.setup(
    name="ecl",
    author="Equinor ASA",
    author_email="fg_sib-scout@equinor.com",
    description="Package for reading and writing the result files from the ECLIPSE reservoir simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/equinor/ecl",
    packages=setuptools.find_packages(
        where="python",
        exclude=["*.tests", "*.tests.*", "tests.*", "tests", "ert.*", "ert"],
    ),
    package_dir={"": "python"},
    license="GPL-3.0",
    platforms="any",
    install_requires=[
        "cwrap",
        "functools32;python_version=='2.7'",
        "future",
        "numpy;python_version>='3.0'",
        "numpy<=1.16.6;python_version=='2.7'",
        "pandas;python_version>='3.0'",
        "pandas<=0.25.3;python_version=='2.7'",
        "six",
    ],
    entry_points={"console_scripts": utility_wrappers()},
    cmake_args=[
        "-DECL_VERSION=" + version,
        "-DBUILD_APPLICATIONS=" + ("ON" if sys.platform == "linux" else "OFF"),
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        # "-DCMAKE_INSTALL_BINDIR=python/ecl/.bin",
        # "-DCMAKE_INSTALL_LIBDIR=python/ecl/.libs",
        # "-DCMAKE_INSTALL_INCLUDEDIR=python/ecl/.include",
        # we can safely pass OSX_DEPLOYMENT_TARGET as it's ignored on
        # everything not OS X. We depend on C++11, which makes our minimum
        # supported OS X release 10.9
        "-DCMAKE_OSX_DEPLOYMENT_TARGET=10.13",
    ],
    # skbuild's test imples develop, which is pretty obnoxious instead, use a
    # manually integrated pytest.
    cmdclass={"test": setuptools.command.test.test},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries",
        "Topic :: Utilities",
    ],
    version=version,
)
