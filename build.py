import ctypes
import os
import platform
import re
import subprocess
import sys
from distutils.version import LooseVersion
from shutil import copyfile, copymode

import distutils.sysconfig as sysconfig

from setuptools import Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name))
        )
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DREQUIRE_OPENMP=0",
            "-DCMAKE_BUILD_TYPE=Debug",
            f"-DPYTHON_INCLUDE_DIR={sysconfig.get_python_inc()}",
            f"-DPYTHON_LIBRARY={sysconfig.get_config_var('LIBDIR')}",
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]
        arch = None
        force_x86 = False

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(
                    cfg.upper(), extdir
                )
            ]
            if (sys.maxsize > 2 ** 32) or (
                8 * ctypes.sizeof(ctypes.c_voidp) == 64
            ):
                arch = "x64"
            else:
                arch = "Win32"

            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        if "MMSEQ_CMAKE_GENERATOR" in os.environ:
            if len(os.environ["MMSEQ_CMAKE_GENERATOR"]) > 0:
                cmake_args += ["-G", os.environ["MMSEQ_CMAKE_GENERATOR"]]
                arch = None

        if "MMSEQ_CMAKE_ARCH" in os.environ:
            if len(os.environ["MMSEQ_CMAKE_ARCH"]) > 0:
                arch = os.environ["MMSEQ_CMAKE_ARCH"]
                if arch == "x86":
                    force_x86 = True

        if force_x86:
            cmake_args += ["-DFORCE_X86=1"]
            build_args += ["-DFORCE_X86=1"]

        if arch is not None:
            cmake_args += ["-A", arch]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print(
            f'The cmake config will use the following args: {" ".join(cmake_args)}'
        )
        print(
            f'The cmake build will use the following args: {" ".join(build_args)}'
        )

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args,
            cwd=self.build_temp,
            env=env,
            stdout=sys.stdout,
            stderr=sys.stdout,
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args,
            cwd=self.build_temp,
            stdout=sys.stdout,
            stderr=sys.stdout,
        )
        # Copy *_test file to tests directory
        # test_bin = os.path.join(self.build_temp, 'unafold_python_test')
        # self.copy_test_file(test_bin)
        print()  # Add an empty line for cleaner output

    def copy_test_file(self, src_file):
        """
        Copy ``src_file`` to ``dest_file`` ensuring parent directory exists.
        By default, message like `creating directory /path/to/package` and
        `copying directory /src/path/to/package -> path/to/package`
        are displayed on standard output. Adapted from scikit-build.
        """
        # Create directory if needed
        dest_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "tests", "bin"
        )
        if dest_dir != "" and not os.path.exists(dest_dir):
            print("creating directory {}".format(dest_dir))
            os.makedirs(dest_dir)

        # Copy file
        dest_file = os.path.join(dest_dir, os.path.basename(src_file))
        print("copying {} -> {}".format(src_file, dest_file))
        copyfile(src_file, dest_file)
        copymode(src_file, dest_file)

# read the contents of README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

datafiles = ["README.md", "CMakeLists.txt", "MANIFEST.in"]
for folder in ["lib", "src", "util", "example", "data", "cmake", "tests"]:
    datafiles = [*datafiles, *[
        os.path.join(d, f)
        for d, folders, files in os.walk(os.path.join(".", folder))
        for f in files
    ]]


def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """

    setup_kwargs.update(
        {
            "long_description": long_description,
            "long_description_content_type": 'text/markdown',
            "ext_modules": [CMakeExtension("mmseqs/mmseqs_native")],
            "cmdclass": {"build_ext": CMakeBuild},
            "zip_safe": False,
            "include_package_data": True,
            "data_files": datafiles,
        }
    )
