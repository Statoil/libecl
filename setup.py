import os
import skbuild
import setuptools
from setuptools_scm import get_version


version = get_version(relative_to=__file__, write_to="python/ecl/version.py")


with open("README.md") as f:
    long_description = f.read()


skbuild.setup(
    name="ecl",
    author="Equinor ASA",
    author_email="fg_sib-scout@equinor.com",
    description="Package for reading and writing the result files from the ECLIPSE reservoir simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/equinor/ecl",
    packages=setuptools.find_packages(where='python', exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_dir={"": "python"},
    license="GPL-3.0",
    platforms="any",
    install_requires=[
        "cwrap",
        "numpy",
        "pandas",
        "six"
    ],
    entry_points={
        "console_scripts": [
            "ecl_pack.x = ecl.bin:ecl_pack",
            "ecl_unpack.x = ecl.bin:ecl_unpack",
        ]
    },
    cmake_args=[
        "-DECL_VERSION=" + version,
        "-DBUILD_APPLICATIONS=ON",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_INSTALL_BINDIR=python/ecl/.bin",
        "-DCMAKE_INSTALL_LIBDIR=python/ecl/.libs",
        "-DCMAKE_INSTALL_INCLUDEDIR=python/ecl/.include",
        # we can safely pass OSX_DEPLOYMENT_TARGET as it's ignored on
        # everything not OS X. We depend on C++11, which makes our minimum
        # supported OS X release 10.9
        "-DCMAKE_OSX_DEPLOYMENT_TARGET=10.9"
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
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries",
        "Topic :: Utilities"
    ],
    version=version
)
