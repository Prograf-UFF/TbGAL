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

Here you find a brief description of the namespaces, macros, classes, functions, procedures, and operators available for the user. All methods are available with C++ and most of them with Python. The detailed documentation is not ready yet.

Contents:

- [Namespaces](#namespaces)
- [Macros](#macros)
- [Classes and Data Types](#classes-and-data-types)
- [Utilities Constants and Functions](#utilities-constants-and-functions)
- [Products and Basic Operations](#products-and-basic-operations)
- [Overloaded Operators](#overloaded-operators)
- [Tools](#tools)
- [Algebra-Specific Declarations](#algebra-specific-declarations)
  - [Signed](#signed)
  - [Euclidean](#euclidean)
  - [Homogeneous/Projective](#homogeneousprojective)
  - [Mikowski/Spacetime](#mikowskispacetime)
  - [Conformal](#conformal)

### Namespaces

In C++, namespaces are declarative regions that provide scope to the names of the types, function, variables, *etc.*, inside it. TbGAL defines the following namespaces.

| Namespace | Description |
| --- | --- |
| `tbgal` | The main namespace that encloses all TbGAL implementations |
| `tbgal::Euclidean1`, `tbgal::Euclidean2`, `tbgal::Euclidean3`, `tbgal::EuclideanD` | The namespace of Euclidean geometric algebra of R<sup>*n*</sup> |
| `tbgal::Homogeneous1`, `tbgal::Homogeneous2`, `tbgal::Homogeneous3`, `tbgal::HomogeneousD` | The namespace of homogeneous/projective geometric algebra of R<sup>*d*</sup> (*n* = *d* + 1) |
| `tbgal::Minkowski1`, `tbgal::Minkowski2`, `tbgal::Minkowski3`, `tbgal::MinkowskiD` | The namespace of Mikowski/spacetime algebra of R<sup>*d*</sup> (*n* = *d* + 2) |
| `tbgal::Conformal1`, `tbgal::Conformal2`, `tbgal::Conformal3`, `tbgal::ConformalD` | The namespace of conformal geometric algebra of R<sup>*d*</sup> (*n* = *d* + 2) |
| `tbgal::SignedPQ` | The namespace of geometric algebras of R<sup>*p, q*</sup> (*n* = *p* + *q*) with metric signatura *(p, q)* |

The `tbgal` namespace also declares a nested `detail` namespace. This is the namespace where the magic happens. Don't touch it!

According to the TbGAL conventions, the root directory for the header files that you will include in your program is the `tbgal` folder. The core operations may be implemented by TbGAL using different libraries, so you have to indicate the one that will be used. So far, Eigen is the only one available. It can be indicated by including the header file `tbgal/using_Eigen.hpp`. Also, the header file for each above-mentioned namespace is its name preceded by `assuming_` and followed by the `.hpp` extension. Putting both conventions together, we have `tbgal/assuming_Euclidean3.hpp`, `tbgal/assuming_Homogeneous3.hpp`, `tbgal/assuming_Minkowski3.hpp`, `tbgal/assuming_Conformal3.hpp`, and so on.

As an example, if you want to use the Eigen-based implementation of TbGAL with conformal geometric algebra of R<sup>3</sup> then you have to put the following instructions among the first lines of your source code:

```cpp
#include <tbgal/using_Eigen.hpp>
#include <tbgal/assuming_Conformal3.hpp>

using namespace tbgal;
using namespace tbgal::Conformal3;
```

In Python, one only have to import the content of the submodule related to the model of geometry:

```python
from tbgal.conformal3 import *
```

### Macros

Optionally, set the following macros before including TbGAL headers in your program to change some conventions of the library.

| Class | Description |
| --- | --- |
| `TBGAL_DEFAULT_SCALAR_TYPE` | Defines the floating-point type assumed as default by the library for scalar values (default is `std::double_t`) |
| `TBGAL_DEFAULT_INDEX_TYPE` | Defines the signed integral type assumed as default by the library for indices (default is `std::int64_t`) |
| `TBGAL_DEFAULT_FLT_TOLERANCE`, `TBGAL_DEFAULT_DBL_TOLERANCE` | Define the tolerances for round-errors while comparing `std::float_t` and `std::double_t` values, respectively |

### Classes and Data Types

The following basic data types are defined in order to assign a meaning to conventional types, like `double`, `int`, and so on.

| Basic Type | Description |
| --- | --- |
| `DefaultScalarType` | The floating point type assumed as default by the library for scalar values (see `TBGAL_DEFAULT_SCALAR_TYPE`) |
| `DefaultIndexType` | The signed integral type assumed as default by the library for indices (see `TBGAL_DEFAULT_INDEX_TYPE`) |

The following classes correspond to the most important structures of TbGAL.

| Class | Description |
| --- | --- |
| `FactoredMultivector<ScalarType, FactoringProductType>` | Concrete class for multivectors enconding a *k*-blade or *k*-versor using the factorization defined by the `FactoringProductType` tag class |
| `GeometricProduct<MetricSpaceType>`, `OuterProduct<MetricSpaceType>` | Tag classes for the `FactoringProductType` |
| `BaseSignedMetricSpace<P, Q [, MaxN]>` | Abstract superclass of classes implementing the `MetricSpaceType` concept |
| `ConformalMetricSpace<D [, MaxD]>`, `EuclideanMetricSpace<N [, MaxN]>`, `HomogeneousMetricSpace<D [, MaxD]>`, `MinkowskiMetricSpace<D [, MaxD]>`, `SignedMetricSpace<P, Q [, MaxN]>` | Concrete classes implementing the `MetricSpaceType` concept |

| Exception Class | Description |
| --- | --- |
| `NotSupportedError` | An exception to report errors related to not implemented features |

The `D`, `P`, `Q`, `MaxD`, and `MaxN` template arguments of the classes implementing the `MetricSpaceType` concept may be set to non-negative integer values in compilation time. As a result, the dimensionality of the vector space will be constant at runtime. The other option is to set them to `Dynamic` if one plans to change the dimensionality of the vector space at runtime. You don't have to worry about that if you are using a model of geometry defined in one of the `tbgal/assuming_[whatever].hpp` headers.

The explicit use of C++ templates while implementing a solution may be overwhelming. For the sake of simplicity, it is strongly recommended to use the `auto` placeholder type specifier (please, refer to the [C++ specification](https://en.cppreference.com/w/cpp/language/auto) for details) whenever possible.

### Utilities Constants and Functions

Here you find some useful functions to assist the implementation of your program.

| Function | Description |
| --- | --- |
| `e(index)` | Returns an unit basis vector |
| `scalar(arg)` | Converts the given numerical value to a scalar factored multivector |
| `vector(coords...)` | Makes a vector with the given set of coordinates |
| `vector(begin, end)` | Makes a vector with the set of coordinates accessed by the iterators |

### Products and Basic Operations

The following tables present a set of basic products and operations from geometric algebra.

| Product | Description |
| --- | --- |
| `dot(arg1, arg2)` | Dot product |
| `gp(arg1, arg2)` | Geometric/Clifford product |
| `hip(arg1, arg2)` | Hestenes inner product |
| `igp(arg1, arg2)` | Inverse geometric/Clifford product (the argument `rhs` must be a versor)  |
| `lcont(arg1, arg2)` | Left contraction |
| `op(arg1, args...)` | Outer/Wedge product |
| `rcont(arg1, arg2)` | Right contraction |
| `sp(arg1, arg2)` | Scalar product |

| Simple Binary Operation | Description |
| --- | --- |
| `addition(arg1, arg2)`, `add(arg1, arg2)` | Addition |
| `subtraction(arg1, arg2)`, `sub(arg1, arg2)` | Subtraction |

| Sign-Change Operation | Description |
| --- | --- |
| `conjugate(arg)` | Clifford conjugation |
| `involute(arg)` | Grade involution |
| `reverse(arg)` | Reversion |
| `unary_minus(arg)` | Unary minus |
| `unary_plus(arg)` | Unary plus |

| Dualization Operation | Description |
| --- | --- |
| `dual(arg)` | Dualization operation |
| `undual(arg)` | Undualization operation |

| Norm-Based Operation | Description |
| --- | --- |
| `rnorm_sqr(arg)` | Squared reverse norm |
| `rnorm(arg)` | Reverse norm |
| `inverse(arg)`, `inv(arg)` | Inverse of the given versor using the squared reverse norm |
| `unit(arg)` | Unit under reverse norm |

| Transformation Operation | Description |
| --- | --- |
| `apply_even_versor(versor, arg)` | Returns the argument transformed by the even versor using the sandwich product |
| `apply_odd_versor(versor, arg)` | Returns the argument transformed by the odd versor using the sandwich product |
| `apply_rotor(rotor, arg)` | Returns the argument transformed by the rotor using the sandwich product |

### Overloaded Operators

TbGAL overload some C++ operators to make the writing of source code closer to the writing of mathematical expressions with geometric algebra.

It is important to notice that the precedence and associativity of C++ operators are different than the one assumed in mathematical functions. For instance, one would expect that the outer/wedge product `^` would be evaluated before the addition operation in the following expression `a + b ^ c`, because product precedes addition in math. However, in C++ the addition operator (`+`) precedes the bitwise XOR operator (`^`), leading to possible mistakes while implementing mathematical procedures (please, refer to the [C++ specification](https://en.cppreference.com/w/cpp/language/operator_precedence) for details). As a result, the resulting expression in this example would be `(a + b) ^ c`. The use of parenthesis is strongly recommended in order to avoid those mistakes. By rewriting the example, `a + (b ^ c)` will guarantee the expected behavior.

| Arithmetic Operator | Description |
| --- | --- |
| `+arg` | Unary plus (same as `unary_plus(arg)`) |
| `-arg` | Unary minus (same as `unary_minus(arg)`) |
| `arg1 + arg2` | Addition (same as `add(arg1, arg2)`) |
| `arg1 - arg2` | Subtraction (same as `sub(arg1, arg2)`) |
| `arg1 * arg2` | Geometric/Clifford product (same as `gp(arg1, arg2)`) |
| `arg1 / arg2` | Inverse geometric/Clifford product (same as `igp(arg1, arg2)`) |
| `arg1 ^ arg2` | Outer/Wedge product (same as `op(arg1, arg2)`) |

| Input/Output Operator | Description |
| --- | --- |
| `os << arg` | Insert formatted output |

### Tools

TbGAL includes a set of useful functions to help developers to write their programs.

| Function | Description |
| --- | --- |
| `default_tolerance<ValueType>()` | Return the standard tolerance value `tol` assumed for the given value type |

| Testing Function | Description |
| --- | --- |
| `is_blade(arg)` | Returns whether the given argument is a blade |
| `is_zero(arg)` | Returns whether the given argument is equal to *zero* |

| Testing Meta-Function | Description |
| --- | --- |
| `is_multivector_v<Type>` | Returns whether the given type is a factored multivector expression |

### Algebra-Specific Declarations

In the following sub-section, you find declarations that are specific of the respective geometric algebra.

#### Signed

Classes and constants of signed geometric algebras of R<sup>*p, q*</sup>. They are available in the following namespace: `tbgal::SignedPQ`.

| Class | Description |
| --- | --- |
| `SignedMetricSpace<P, Q [, MaxN]>` | Orthogonal metric space with signature (*p*, *q*) (*n* = *p* + *q*) |

| Constant Value | Description |
| --- | --- |
| `SPACE` | An instance of the orthogonal metric space class with signature (*p*, *q*) |

#### Euclidean

Classes, constants, functions, and operations of Euclidean geometric algebras of R<sup>*n*</sup>. They are available in the following namespaces: `tbgal::Euclidean1`, `tbgal::Euclidean2`, `tbgal::Euclidean3`, and `tbgal::EuclideanD`.

| Class | Description |
| --- | --- |
| `EuclideanMetricSpace<N [, MaxN]>` | Euclidean metric space |

| Constant Value | Description |
| --- | --- |
| `e1`, `e2`, ..., `eN` | Euclidean basis vector (same as `e(1)`,  `e(2)`, ..., `e(N)`) |
| `SPACE` | An instance of the Euclidean metric space class |

| Function | Description |
| --- | --- |
| `euclidean_vector(coords...)` | Makes an Euclidean vector with the given set of coordinates |
| `euclidean_vector(begin, end)` | Makes an Euclidean vector with the set of coordinates accessed by the iterators |

#### Homogeneous/Projective

Classes, constants, functions, and operations of homogeneous/projective geometric algebras of R<sup>*d*</sup> (*n* = *d* + 1). They are available in the following namespaces: `tbgal::Homogeneous1`, `tbgal::Homogeneous2`, `tbgal::Homogeneous3`, and `tbgal::HomogeneousD`.

| Class | Description |
| --- | --- |
| `HomogeneousMetricSpace<D [, MaxD]>` | Homogeneous/Projective metric space |

| Constant Value | Description |
| --- | --- |
| `e1`, `e2`, ..., `eD` | Euclidean basis vector (same as `e(1)`,  `e(2)`, ..., `e(D)`) |
| `ep` | Positive extra basis vector interpreted as the point at the origin (same as `e(D + 1)`) |
| `SPACE` | An instance of the homogeneous/projective metric space class |

| Function | Description |
| --- | --- |
| `direction(coords...)` | Makes a direction vector using the given set of coordinates |
| `direction(begin, end)` | Makes a direction vector using the set of coordinates accesses by the iterators |
| `euclidean_vector(coords...)` | Makes an Euclidean vector with the given set of coordinates |
| `euclidean_vector(begin, end)` | Makes an Euclidean vector with the set of coordinates accessed by the iterators |
| `point(coords...)` | Makes an unit point using the given set of coordinates |
| `point(begin, end)` | Makes an unit point using the set of coordinates accesses by the iterators |

#### Mikowski/Spacetime

Classes, constants, functions, and operations of Mikowski/spacetime geometric algebras of R<sup>*d*</sup> (*n* = *d* + 2). They are available in the following namespaces: `tbgal::Minkowski1`, `tbgal::Minkowski2`, `tbgal::Minkowski3`, and `tbgal::MinkowskiD`.

| Class | Description |
| --- | --- |
| `MinkowskiMetricSpace<D [, MaxD]>` | Minkowski/Spacetime metric space |

| Constant Value | Description |
| --- | --- |
| `e1`, `e2`, ..., `eD` | Euclidean basis vector (same as `e(1)`, `e(2)`, ..., `e(D)`) |
| `ep` | Positive extra basis vector (same as `e(D + 1)`) |
| `em` | Negative extra basis vector (same as `e(D + 2)`) |
| `SPACE` | An instance of the Minkowski/spacetime metric space class |

| Function | Description |
| --- | --- |
| `euclidean_vector(coords...)` | Makes an Euclidean vector with the given set of coordinates |
| `euclidean_vector(begin, end)` | Makes an Euclidean vector with the set of coordinates accessed by the iterators |
| `point(coords...)` | Makes an unit point using the given set of coordinates |
| `point(begin, end)` | Makes an unit point using the set of coordinates accesses by the iterators |

#### Conformal

Classes, constants, functions, and operations of conformal geometric algebras of R<sup>*d*</sup> (*n* = *d* + 2). They are available in the following namespaces: `tbgal::Conformal1`, `tbgal::Conformal2`, `tbgal::Conformal3`, and `tbgal::ConformalD`.

| Class | Description |
| --- | --- |
| `ConformalMetricSpace<D [, MaxD]>` | Conformal metric space |

| Constant Value | Description |
| --- | --- |
| `e1`, `e2`, ..., `eD` | Euclidean basis vector (same as `e(1)`,  `e(2)`, ..., `e(D)`) |
| `no` | Null vector interpreted as the point at the origin (same as `e(D + 1)`) |
| `ni` | Null vector interpreted as the point at infinity (same as `e(D + 2)`) |
| `SPACE` | An instance of the conformal metric space class |

| Function | Description |
| --- | --- |
| `euclidean_vector(coords...)` | Makes an Euclidean vector with the given set of coordinates |
| `euclidean_vector(begin, end)` | Makes an Euclidean vector with the set of coordinates accessed by the iterators |
| `point(coords...)` | Makes an unit point using the given set of coordinates |
| `point(begin, end)` | Makes an unit point using the set of coordinates accesses by the iterators |

## 6. Related Project

Please, visit the GitHub repository of the [**ga-benchmark**](https://github.com/ga-developers/ga-benchmark) project for a benchmark comparing the most popular libraries, library generators, and code optimizers for geometric algebra.

## 7. License

This software is licensed under the GNU General Public License v3.0. See the [`LICENSE`](LICENSE) file for details.
