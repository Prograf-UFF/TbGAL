# TbGAL: Tensor-Based Geometric Algebra Library

TbGAL is a C++/Python library for Euclidean, homogeneous/projective, Mikowski/spacetime, conformal, and arbitrary [geometric algebras](https://en.wikipedia.org/wiki/Geometric_algebra) with assuming (p, q) metric signatures.

Geometric algebra is a powerful mathematical system encompassing many mathematical concepts (*e.g.*, [complex numbers](https://en.wikipedia.org/wiki/Complex_number), [quaternions algebra](https://en.wikipedia.org/wiki/Quaternion_algebra), [Grassmann-Cayley algebra](https://en.wikipedia.org/wiki/Grassmann%E2%80%93Cayley_algebra), and [Plücker coordinates](https://en.wikipedia.org/wiki/Pl%C3%BCcker_coordinates)) under the same framework. Geometric algebra is mainly based on the algebraic system called [Clifford algebra](https://en.wikipedia.org/wiki/Clifford_algebra), but with a strong emphasis on geometric interpretation. In geometric algebra, subspaces are treated as primitives for computation. As such, it is an appropriate mathematical tool for modeling and solving geometric problems in physics, chemistry, engineering, and computer science.

TbGAL represents blades (and versors) in their decomposed state as the outer product (and geometric product), rather than using their representation as a weighted summation of basis blades in ⋀ℝ<sup>*n*</sup>. The main advantage of the factorized approach is that it is able to compute geometric algebra operations in higher dimensions, *i.e.*, assume multivectors space ⋀ℝ<sup>*n*</sup> with *n* > 256. In terms of memory, TbGAL stores only 1 + *n*<sup>2</sup> coefficients per blade or versor in worst case, while operations have maximum complexity of O(*n*<sup>3</sup>).

Please cite our [Advances in Applied Clifford Algebras](https://doi.org/10.1007/s00006-020-1053-1) paper if you use this code in your research. The paper presents a complete description of the library:

```{.bib}
@Article{sousa_fernandes-aaca-30(2)-2020,
  author  = {Sousa, Eduardo V. and Fernandes, Leandro A. F.},
  title   = {{TbGAL}: a tensor-based library for geometric algebra},
  journal = {Advances in Applied Clifford Algebras},
  year    = {2020},
  volume  = {30},
  number  = {2},
  pages   = {75},
  doi     = {https://doi.org/10.1007/s00006-020-1053-1},
  url     = {https://github.com/Prograf-UFF/TbGAL},
}
```

Please, let Eduardo Vera Sousa ([http://www.ic.uff.br/~eduardovera](http://www.ic.uff.br/~eduardovera)) and Leandro A. F. Fernandes ([http://www.ic.uff.br/~laffernandes](http://www.ic.uff.br/~laffernandes)) know if you want to contribute to this project. Also, do not hesitate to contact them if you encounter any problems.

**Contents:**

1. [Requirements](#1-requirements)
2. [How to Build and Install TbGAL](#2-how-to-build-and-install-tbgal)
3. [Compiling and Running Examples](#3-compiling-and-running-examples)
4. [Compiling and Running Unit Tests](#4-compiling-and-running-unit-tests)
5. [Documentation](#5-documentation)
6. [Related Project](#6-related-project)
7. [License](#7-license)

## 1. Requirements

Make sure that you have the following tools before attempting to use TbGAL.

Required tools:

- Your favorite [C++17](https://en.wikipedia.org/wiki/C%2B%2B17) compiler.
- [CMake](https://cmake.org) (version >= 3.14) to automate installation and to build and run examples.

Required C++ library:

- [Eigen](http://eigen.tuxfamily.org) (version >= 3) to evaluate basic matrix algebra routines.

Required tools, if you want to use TbGAL with Python:

- [Python 2 or 3](https://www.python.org) interpreter, if you want to build and use TbGAL with Python.

Required Python packages and C++ libraries, if you want to use TbGAL with Python:

- [NumPy](https://numpy.org), the fundamental package for scientific computing with Python.
- [Boost.Python](https://www.boost.org/doc/libs/release/libs/python/doc/html/index.html) (version >= 1.56), a C++ library which enables seamless interoperability between C++ and the Python programming language.
- [Boost.NumPy](https://www.boost.org/doc/libs/release/libs/python/doc/html/numpy/index.html), a C++ library that extends Boost.Python to NumPy.

Required C++ libraries, if you want to run unit tests:

- [GATL](https://github.com/laffernandes/gatl) is another C++ library for geometric algebra.
- [Google Test](https://github.com/google/googletest) is a unit testing library for the C++ programming language, based on the xUnit architecture.

Optional tool to use TbGAL with Python:

- [Virtual enviroment](https://wiki.archlinux.org/index.php/Python/Virtual_environment) to create an isolated workspace for a Python application.

## 2. How to Build and Install TbGAL

Use the [git clone](https://git-scm.com/docs/git-clone) command to download the project, where `<tbgal-dir>` must be replaced by the directory in which you want to place TbGAL's source code, or remove `<tbgal-dir>` from the command line to download the project to the `./TbGAL` directory:

```bash
git clone https://github.com/Prograf-UFF/TbGAL.git <tbgal-dir>
```

TbGAL/C++ is a pure template library. Therefore, there is no binary library to link to, but you have to install the header files. Optionally, you can build and install the TbGAL/Python front-end.

Use CMake to copy TbGAL's header files to the common include directory in your system (*e.g.*, `/usr/local/include`, in Linux). The basic steps for installing TbGAL/C++ using CMake look like this in Linux:

```bash
cd <tbgal-dir>
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
```

Notice that you may use the `-G <generator-name>` option of CMake's command-line tool to choose the build system (*e.g.*, Unix makefiles, Visual Studio, etc.). Please, refer to [CMake's Help](https://cmake.org/cmake/help/latest/manual/cmake.1.html) for a complete description of how to use the CMake's command-line tool.

After installation, CMake will find TbGAL/C++ using the command `find_package(TbGAL)` (see [CMake documentation](https://cmake.org/cmake/help/latest/command/find_package.html) for details). In addition, you will be able to use the `TbGAL_INCLUDE_DIRS` variable in the `CMakeList.txt` file of your program while defining the include directories of your C++ project or targets.

TbGAL/Python is a back-end to access TbGAL/C++ from a Python environment. In this case, you have to build and install the TbGAL/Python modules using the commands presented bellow:

```bash
cmake --build . --config Release --parallel 8 --target install
```

It is important to emphasize that both Python 2 and 3 are supported. Please, refer to [CMake's documentation](https://cmake.org/cmake/help/latest/module/FindPython.html) for details about how CMake finds the Python interpreter, compiler, and development environment.

Finally, add `<cmake-install-prefix>/lib/tbgal/python/<python-version>` to the the `PYTHONPATH` environment variable. The `<cmake-install-prefix>` placeholder usually is `/usr/local` on Linux, and `C:/Program Files/TbGAL` or `C:/Program Files (x86)/TbGAL` on Windows. But it may change according to what was set in CMake. The `<python-version>` placeholder is the version of the Python interpreter found by CMake.

Set the `PYTHONPATH` variable by calling following command in Linux:

```bash
export PYTHONPATH="$PYTHONPATH:<cmake-install-prefix>/lib/tbgal/python/<python-version>"
```

But this action is not permanent. The new value of `PYTHONPATH` will be lost as soon as you close the terminal. A possible solution to make an environment variable persistent for a user's environment is to export the variable from the user's profile script:

  1. Open the current user's profile (the `~/.bash_profile` file) into a text editor.
  2. Add the export command for the `PYTHONPATH` environment variable at the end of this file.
  3. Save your changes.

Execute the following steps to set the `PYTHONPATH` in Windows:

  1. From the *Windows Explorer*, right click the *Computer* icon.
  2. Choose *Properties* from the context menu.
  3. Click the *Advanced system settings* link.
  4. Click *Environment Variables*. In the section *System Variables*, find the `PYTHONPATH` environment variable and select it. Click *Edit*. If the `PYTHONPATH` environment variable does not exist, click *New*.
  5. In the *Edit System Variable* (or *New System Variable*) window, specify the value of the `PYTHONPATH` environment variable to include `"<cmake-install-prefix>/lib/tbgal/python/<python-version>"`. Click *OK*. Close all remaining windows by clicking *OK*.
  6. Reopen yout Python environment.

## 3. Compiling and Running Examples

The basic steps for configuring and building the C++ example of the TbGAL look like this:

```bash
cd <tbgal-dir>/cpp/example
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release --parallel 8
```

Call the executables files placed at `<tbgal-dir>/cpp/example/build` on Linux and at `<tbgal-dir>\cpp\example\build\Release` on Windows.

Recall that `<tbgal-dir>` is the directory in which you placed TbGAL's source code.

Use the files in the `<tbgal-dir>/cpp/example` directory as examples of how to use TbGAL in your C++ program. For instance, after installation of the TbGAL library, CMake will find TbGAL using the command `find_package(TbGAL)` (see the [`<tbgal-dir>/cpp/example/CMakeLists.txt`](cpp/example/CMakeLists.txt) file and the [CMake documentation](https://cmake.org/cmake/help/latest/command/find_package.html) for details). Also, you will be able to use the `TbGAL_INCLUDE_DIRS` variable in the `CMakeList.txt` file of your program while defining the include directories of your C++ project or targets. In your source code, you have to use the `#include <tbgal/using_Eigen.hpp>` directive to instrument the library to perform matrix computations using [Eigen](http://eigen.tuxfamily.org) and the `#include <tbgal/assuming_[some-model].hpp>` directive to assume some pre-defined model of geometry.

Similarly, you will find examples of how to use the TbGAL library with Python in the [`<tbgal-dir>/python/example/py2`](python/example/py2) and [`<tbgal-dir>/python/example/py3`](python/example/py3) directories.

## 4. Compiling and Running Unit Tests

The basic steps for configuring and building the C++ unit tests of the TbGAL look like this:

```bash
cd <tbgal-dir>/cpp/test
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --parallel 8
```

The configuration process will donwload, build, and install [Google Test](https://github.com/google/googletest) automatically. [GATL](https://github.com/laffernandes/gatl) must be installed by you first. If successfull, call the `test` target to run all tests:

```bash
cmake --build . --target test
```

## 5. Documentation

Documentation is under construction. Sorry!

## 6. Related Project

Please, visit the GitHub repository of the [**ga-benchmark**](https://github.com/ga-developers/ga-benchmark) project for a benchmark comparing the most popular libraries, library generators, and code optimizers for geometric algebra.

## 7. License

This software is licensed under the GNU General Public License v3.0. See the [`LICENSE`](LICENSE) file for details.
