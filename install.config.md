> See [install.md](./install.md) for general installation guidelines.

This file contains information about the ParaMonte library build configuration flags.

The ParaMonte library build scripts accept a large number of
optional flags that can be used to configure the library build.

Nearly all configuration flags can be specified as command line arguments
to the [install.bat](./install.bat.md), [install.sh](./install.sh.md), or
the [CMake binary executable](./CMakeLists.md). Some flags may be essential
to directly call the CMake executable binary, while others are available
only as optional arguments to the library install scripts mentioned above.

Due to the sheer number of available optional flags, they are categorized
in tiers sorted by their relevance and importance to the end users.

1.  [TIER-1 ParaMonte library build configuration flags](#TIER-1-ParaMonte-library-build-configuration-flags)
2.  [TIER-2 ParaMonte library build configuration flags](#TIER-2-ParaMonte-library-build-configuration-flags)
3.  [TIER-3 ParaMonte library build configuration flags](#TIER-3-ParaMonte-library-build-configuration-flags)
4.  [TIER-4 ParaMonte library build configuration flags](#TIER-4-ParaMonte-library-build-configuration-flags)

## TIER-1 ParaMonte library build configuration flags

The ParaMonte TIER-1 build flags set the most important and useful build configurations.
Assuming the required compilers and external (e.g., Coarray, MPI parallel) libraries are
available on the system, these flags can readily customized the library build.

### `build`

Specifies the library build type for which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --build "build_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dbuild="build_type"
    ```

where `build_type` can be any of the following values.

Value               | Usage
--------------------|------
`debug`             | For debugging.
`testing`           | For quick testing with minimal optimizations.
`release`           | For production-ready application builds (default).
`ipo`               | Same as specifying `release` but with link-time optimizations.
`tuned`             | Same as specifying `ipo` but builds for all supported architectures.
`native`            | Same as specifying `ipo` but builds only for the current architecture (**non-portable**).

> When specified for install scripts, the `build` flag can take a semi-colon-separated
list of library build values for which the library will be built in the specified order.

> If `testing` is specified, it is automatically converted to `RelwithDebInfo`.

> This option, if provided, will overwrite the value of `CMAKE_BUILD_TYPE` in CMake scripts.

> Beware that `ipo`, `tuned`, `native` builds are currently only recognized
used with the `intel` and `gnu` Fortran compilers. If specified with unsupported
compilers, they are automatically converted to `release`.

> The `ipo`, `tuned`, and `native` builds can take very long and generate large binaries.

**optional**. The default value for `build_type` is `release`.

### `lang`

Specifies the target programming language(s) for
which the library will be built and accessed from.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --lang "programming_language"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dlang="programming_language"
    ```

where `programming_language` can be any of the following
values (which defines the corresponding macro in the library),

Value               | Target programming language       | Library macro definition
--------------------|-----------------------------------|-------------------------
`c`                 | C                                 | `C_ENABLED=1`
`cpp`               | C++                               | `CPP_ENABLED=1`
`fortran`           | Fortran                           | `FORTRAN_ENABLED=1`
`go`                | Go                                | `GO_ENABLED=1`
`java`              | Java                              | `JAVA_ENABLED=1`
`julia`             | Julia                             | `JULIA_ENABLED=1`
`mathematica`       | Mathematica                       | `MATTHEMATICA_ENABLED=1`
`matlab`            | MATLAB                            | `MATLAB_ENABLED=1`
`python`            | Python                            | `PYTHON_ENABLED=1`
`r`                 | R                                 | `R_ENABLED=1`

> When specified for install scripts, the `lang` flag can take a semi-colon-separated
list of programming languages for which the library will be built in the specified order.

**optional**. The default value for `programming_language` is `fortran`.

### `lib`

Specifies the library file type for which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --lib "library_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dlib="library_type"
    ```

where `library_type` can be any of the following values.

Value               | Usage
--------------------|------
`static`            | For static library file (`*.lib`, `*.a` generation) (and linking).
`shared`            | For dynamic library file (`*.dll, *.dylib`, `*.so`) generation (and linking).

> When specified for install scripts, the `lib` flag can take a semi-colon-separated
list of possible values for which the library will be built in the specified order.

**optional**. The default value for `library_type` is `shared`.

### `mem`

Specifies the library memory usage type for which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --mem "memory_allocation_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dmem="memory_allocation_type"
    ```

where `memory_allocation_type` can be any of the following values.

Value               | Usage
--------------------|------
`stack`             | All arrays (, if possible, including allocatable,) are stored on the Stack
`heap`              | Most arrays of significant size (e.g., > 10 elements) will be stored on the Heap.

> The use of Heap memory can be significantly slower than the Stack in isolated cases.

> The stack memory allocations can lead to "stack overflow" or "segmentation fault"
that are frequently extremely hard to identify.

> If you use `stack` be sure to increase the stack memory allocated
by the Operating System for your output binary before execution.

> When specified for install scripts, the `mem` flag can take a semi-colon-separated
list of possible values for which the library will be built in the specified order.

**optional**. The default value for `memory_allocation_type` is `heap`.

### `par`

Specifies the library parallelization paradigm for which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --par "parallelization_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dpar="parallelization_type"
    ```

where `parallelization_type` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | No parallelism enabled (serial).
`serial`            | No parallelism enabled (serial).
`mpi`               | MPI parallelism enabled.
`omp`               | OpenMP parallelism enabled.
`openmp`            | OpenMP parallelism enabled.
`cafsingle`         | Coarray parallelism syntax enabled (but effectively serial).
`cafshared`         | Coarray parallelism enabled on shared memory architecture.
`cafdist`           | Coarray parallelism enabled on distributed memory architecture.

> When specified for install scripts, the `par` flag can take a semi-colon-separated
list of possible values for which the library will be built in the specified order.

> Beware that the Coarray parallelism in the current library version 2.0 is non-functional.

**optional**. The default value for `parallelization_type` is `none`.

## TIER-2 ParaMonte library build configuration flags

The ParaMonte TIER-2 build flags set the additional optional arguments that
run library benchmarks, examples, or tests, or are critical for the correct
selection of compilers or external compilation libraries, or customize the
build and installation folders.

### `bdir`

Specifies the library build directory where all
CMake and other relevant files will be stored
before outputting the final product to the
specified deployment directory via `ddir`.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --bdir "cmake_build_directory_path"
    ```
+   Usage (with `cmake` binary executable): **Not Available.**

where `cmake_build_directory_path` is the path
to the directory where the library should be built.

> This option cannot be specified as a flag to the CMake binary
> because the CMake build scripts must already exist in this folder.

> **WARNING**
> The specified value for this flag must be an **absolute** (full) path.
> Relative paths currently result in library build failure.

**optional**. The default value for `cmake_build_directory_path` is,

```bash
bdir="${root}/bld/${os}/${arch}/${csid}${csvs}/${build}/${lib}/${mem}/${par}/${checking}/${lang}"
```

where

+   `${root}`   is replaced with the path to the root directory of the project where the LICENSE file exists.
+   `${os}`     is replaced with the lower-case operating system name (`windows`, `linux`, `darwin`, `mingw`, `msys`, ...),
+   `${arch}`   is replaced with the lower-case CPU architecture (`amd64`, `arm64`, ...),
+   `${csid}`   is replaced with the lower-case compiler vendor name (`intel`, `gnu`, ...),
+   `${csvs}`   is replaced with the lower-case compiler vendor major version,
+   `${bld}`    is replaced with the specified value for the `build` configuration flag (`debug`, `release`, ...),
+   `${lib}`    is replaced with the specified value for the `lib` configuration flag (`static`, `shared`, ...),
+   `${mem}`    is replaced with the specified value for the `mem` configuration flag (`stack`, `heap`, ...),
+   `${par}`    is replaced with a value determined from the `mem` configuration flag:

    Value               | Scenario
    --------------------|--------------------------------------------------------------
    `cafsingle`         | If the library is built for Coarray single-image parallelism.
    `cafshared`         | If the library is built for Coarray shared-memory parallelism.
    `cafdist`           | If the library is built for Coarray distributed-memory parallelism.
    `mpi`               | If the library is built for MPI parallelism using an unknown MPI distribution.
    `impi`              | If the library is built for MPI parallelism using an Intel MPI distribution.
    `mpich`             | If the library is built for MPI parallelism using an MPICH MPI distribution.
    `openmpi`           | If the library is built for MPI parallelism using an OpenMPI distribution.
    `openmp`            | If the library is built for OpenMP parallelism.
    `serial`            | If the library is built for serial applications.

### `bench`

Specifies the library benchmarks to build and run after building library.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --bench "benchmark_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dbench="benchmark_list"
    ```

where `benchmark_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | No library benchmarks will be built.
`all`               | Build and run all library benchmarks.
file                | Build and run all library benchmarks listed in the file.
b1;b2;...;bn        | Build and run all modules and procedures benchmarks matching the semicolon-separated strings.

> The specified file can be automatically generated by the `getList.py`
> Python script in the corresponding language subfolder of the `benchmark`
> folder in the root directory of the project.

> Specified items in a given file can be excluded by
> adding an exclamation mark at the beginning of the line.

> Setting this option to anything other than `none` will lead to building
> and running the corresponding ParaMonte library benchmarks on the command line.

> The benchmarks will be built and run in the corresponding library benchmark subfolder in the build folder.

> Not all supported programming languages may have benchmarks.
> Where unavailable, this option should not be specified.

**optional**. The default value for `benchmark_list` is `none`.

### `blas`

Specifies the BLAS implementation against which certain ParaMonte library routines will be linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --blas "cmake_blas_vendor"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dblas="cmake_blas_vendor"
    ```

where `cmake_blas_vendor` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | No BLAS library will be linked to the ParaMonte library.
`default`           | The local OpenBLAS submodule within the ParaMonte repository will be built and used for linking.
`any`               | Any BLAS library CMake finds will be linked to the ParaMonte library.
vendor              | Any vendor name recognized by the CMake `BLA_VENDOR`.

> When the value `default` is specified, the BLAS library will be built
> via the default OpenBLAS local external submodule of the ParaMonte library,
> if the OpenBLAS submodule exists in the `external` folder of the local clone
> of the ParaMonte repository. If the local copy of `OpenBLAS` library exists,
> the library will be linked against OpenBLAS, otherwise, the library will be
> built without any BLAS dependency.

> This option relevant to the ParaMonte library benchmarks and linear algebra routines.

**optional**. The default value for `cmake_blas_vendor` is `none`.

### `checking`

Specifies the library runtime checking policy for which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --checking "checking_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dchecking="checking_type"
    ```

where `checking_type` can be any of the following values.

Value               | Usage
--------------------|------
`nocheck`           | No runtime checks will be performed.
`checked`           | All runtime checks will be performed.

> When specified for install scripts, the `checking` flag can take a semi-colon-separated
list of possible values for which the library will be built in the specified order.

> Specifying `checked` leads to defining the FPP macro `CHECK_ENABLED=1`,
activating additional diagnostic runtime checks within the library.

> This option is very useful for debugging purposes but significantly
degrades the runtime performance and increases the library size.

**optional**. The default value for `checking_type` is `nocheck`.

### `ddir`

Specifies the library deployment directory to
which the full ParaMonte package will be copied.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --ddir "deploy_directory_path"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dddir="deploy_directory_path"
    ```

where `deploy_directory_path` is the path to the
directory where the library binaries should be deployed.

**optional**. The default value for `deploy_directory_path` is
`"./bin"` where `.` refers to the path to the root directory
of the project where the main `CMakeLists.txt` file exists.

### `exam`

Specifies the library examples to build and run after building library.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --exam "example_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dexam="example_list"
    ```

where `example_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | No library examples will be built.
`all`               | Build and run all library examples.
file                | Build and run all library examples listed in the file.
e1;e2;...;en        | Build and run all modules and procedures examples matching the semicolon-separated strings.

> The specified file can be automatically generated by the `getList.py`
> Python script in the corresponding language subfolder of the `example`
> folder in the root directory of the project.

> The file path can be specified either as an absolute path or a path relative
> to the root directory of the ParaMonte repository on the local machine.

> Specific items in a given file can be excluded by
> adding an exclamation mark at the beginning of the line.

> Setting this option to anything other than `none` will lead to building
> and running the corresponding ParaMonte library examples on the command line.

> The examples will be built and run in the corresponding library example subfolder in the build folder.

**optional**. The default value for `example_list` is `none`.

### `fc`

Specifies the path to the Fortran compiler binary
executable file with which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --fc "fortran_compiler_executable_path"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dfc="fortran_compiler_executable_path"
    ```

where `fortran_compiler_executable_path` points to the
specific Fortran compiler executable name (or its path) to be used.

> Setting this will overwrite the default CMake compiler choice by setting
`CMAKE_Fortran_COMPILER` CMake variable **before** project initiation.

> We highly recommend to specify the Fortran compiler choice explicitly via this argument
as CMake often has difficultly choosing the right Fortran compiler among several options.

**optional**. The default value for `fortran_compiler_executable_path`
is automatically determined by CMake or the build scripts.

### `lapack`

Specifies the LAPACK implementation against which certain ParaMonte library routines will be linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --lapack "cmake_lapack_vendor"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dlapack="cmake_lapack_vendor"
    ```

where `cmake_lapack_vendor` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | No LAPACK library will be linked to the ParaMonte library.
`default`           | The local OpenBLAS submodule within the ParaMonte repository will be built and used for linking.
`any`               | Any LAPACK library CMake finds will be linked to the ParaMonte library.
vendor              | Any vendor name recognized by the CMake `BLA_VENDOR`.

> When the value `default` is specified, the LAPACK library will be built
> via the default OpenBLAS local external submodule of the ParaMonte library,
> if the OpenBLAS submodule exists in the `external` folder of the local clone
> of the ParaMonte repository. If the local copy of `OpenBLAS` library exists,
> the library will be linked against OpenBLAS, otherwise, the library will be
> built without any LAPACK dependency.

> This option relevant to the ParaMonte library benchmarks and linear algebra routines.

**optional**. The default value for `cmake_lapack_vendor` is `none`.

### `matlabdir`

Specifies the path to the **root directory** of MATLAB against which the library will be linked.
Within a MATLAB session, this directory path is returned by the MATLAB intrinsic function `matlabroot`.
An example such path returned in MATLAB 2024a environment is: `'C:\Program Files\MATLAB\R2024a'` on Windows platforms.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --matlabdir "matlab_root_dir_path"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dmatlabdir="matlab_root_dir_path"
    ```

where `matlab_root_dir_path` is the path to the root directory of MATLAB.

> This option is relevant only to builds where the `lang` configuration flag is set to `matlab`.
> The specified value is ignored for all language builds of the ParaMonte library.

> Specifying this option helps CMake find the MATLAB library dependencies.

> If specified, it sets the value of the CMake configuration variable `Matlab_ROOT_DIR`.

> This option is essential when the library is to be linked against
> a particular installation of MATLAB among multiple installations.

> We highly recommend to specify the MATLAB choice explicitly via this argument
as CMake often has difficultly choosing the right MATLAB version among several options.

> **NOTE**
> If you are a ParaMonte developer and aim to specify a MATLAB installation against which you intend to build the library,
> always ensure to install the oldest possible compatible MATLAB version. This ensures the generated MEX files are compatible
> with all newer MATLAB version. The opposite does not generally hold. For example, MATLAB 2024a MEX files are not compatible with 
> MATLAB 2020b environment and libraries. 

**optional**. The default value for `matlab_root_dir_path` is set automatically by CMake.

### `me`

Specifies the path to the MPI launcher `mpiexec`
binary executable file for MPI-parallel applications.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --me "mpiexec_path"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dme="mpiexec_path"
    ```

where `mpiexec_path` is the path (or name) of the `mpiexec` executable to be used.

> This option is relevant only to MPI-parallel library builds and is otherwise ignored.

> The specified value will the overwrite the CMake MPI library choice by
setting `MPIEXEC_EXECUTABLE` CMake variable **before** project initiation.

> We highly recommend to specify the MPI library choice explicitly via this argument
as CMake often has difficultly choosing the right MPI compiler among several options.

**optional**. The default value for `mpiexec_path`
is automatically determined by CMake or the build scripts.

### `test`

Specifies the library testing mode for which the library will be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --test "testing_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dtest="testing_type"
    ```

where `testing_type` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | The build is configured for no testing.
`all`               | The build activates all library tests.

> **WARNING**
> Enabling tests requires setting the value of the build flag `mod` to `"all"`.
> This is particularly important for the programming languages other than Fortran
> where the default value of `mod` is usually not `"all"`.

> **WARNING**
> The ParaMonte extended precision tests are prone to failure.
> This is due to GNU compiler bugs for extended precision arithmetic.
> To avoid bug-induced test failures when using GNU compilers,
> you can additionally specify the [`--rki "1;2"`](#rki) build
> to build the library and its test for only the
> single and double `real` type precisions.

**optional**. The default value for `testing_type` is `none`.

## TIER-3 ParaMonte library build configuration flags

The ParaMonte TIER-3 build flags are mostly relevant to the ParaMonte library developers
or advanced users who wish to further customize the library configurations for their needs.
Some of the available options with intricate implications require careful attention before usage.

### `benchpp`

Specifies the library benchmarks postprocessing scripts to run after building the library.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --benchpp "benchmark_postprocessing_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dbenchpp="benchmark_postprocessing_list"
    ```

where `benchmark_postprocessing_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | Do not run benchmarks postprocessing scripts.
`all`               | Run all benchmarks postprocessing scripts.
file                | Run all library benchmarks postprocessing scripts listed in the file.
b1;b2;...;bn        | Run all modules and procedures benchmarks postprocessing scripts matching the semicolon-separated strings.

> The specified file can be automatically generated by the `getList.py`
> Python script in the corresponding language subfolder of the `benchmark`
> folder in the root directory of the project.

> Specific items in a given file can be excluded by
> adding an exclamation mark at the beginning of the line.

> Setting this option to anything other than `none` will lead to running
the Python postprocessing scripts for the corresponding ParaMonte library
benchmarks that have been previously built and run on the command line.

> Frequently, postprocessing scripts merely generate visualizations.
> The visualizations are essential for the ParaMonte library documentation.

> Setting this option to anything other than `none` will lead to running
> the corresponding ParaMonte library benchmarks postprocessing scripts on the command line.

> The benchmarks postprocessing scripts will be executed in the corresponding library benchmark subfolder in the build folder.

> Not all supported programming languages may have benchmarks postprocessing scripts.
> Check the `benchmark` folder in the root directory of the repository before using this flag.

**optional**. The default value for `benchmark_postprocessing_list`
is the same as the value set for the option `bench`.

### `cfi`

Specifies whether the library must be built with C-Fortran interoperable types and kinds.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --cfi "c_fortran_interoperability_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dcfi="c_fortran_interoperability_list"
    ```

where `c_fortran_interoperability_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | Sets the macro `CFI_ENABLED=0` to disable C-Fortran type interoperability.
`all`               | Sets the macro `CFI_ENABLED=1` to enable C-Fortran type interoperability.

> This is a low level option and mostly useful to the developers of the ParaMonte library.

> This option has extremely limited usage for the end users and may be removed in the future.

> The `none` value can be specified only when `lang` is set to `fortran`.

**optional**. The default is `none` if `lang` is set to `fortran`, otherwise `all`.

### `codecov`

Determines whether Code Coverage report must be generated.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --codecov "code_coverage_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dcodecov="code_coverage_type"
    ```

where `code_coverage_type` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | Sets the macro `CODECOV_ENABLED=0` to disable Code Coverage analysis.
`all`               | Sets the macro `CODECOV_ENABLED=1` to enable Code Coverage analysis.

> This is a low level option and mostly useful to the developers of the ParaMonte library.

> This option is currently only tested and verified to work on Linux systems.
> This option is currently only tested and verified to work with GNU compilers.
> This option is currently only tested and verified to work with GNU GCOV and LCOV software.
> This option is currently only tested and verified to work with the `install.sh` build script.

> **WARNING**
> Generating Code Coverage report currently requires the GNU Fortran compiler `>10.3` and a **compatible-version** of **gcov** and **lcov** software.
> The coverage report generation is bound to fail if the compiler and GCOV software version are incompatible.

**optional**. The default value for `codecov` is `none`.

### `deps`

Specifies the library dependencies to copy to the final deployment directory.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --deps "dependencies_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Ddeps="dependencies_list"
    ```

where `dependencies_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | The lower-level library dependencies will be NOT copied to the installation location.
`all`               | The lower-level library dependencies will be copied to the installation location.

> This is a low level option and mostly useful to the developers of the ParaMonte library.

> **NOTE** This flag is currently non-functional on Windows platforms possibly due to a CMake bug related to `TARGET_SONAME_FILE`.

**optional**. The default value for `dependencies_list` is `none`.

### `exampp`

Specifies the library examples postprocessing scripts to run after building the library.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --exampp "example_postprocessing_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dexampp="example_postprocessing_list"
    ```

where `example_postprocessing_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | Do not run examples postprocessing scripts.
`all`               | Run all examples postprocessing scripts.
file                | Run all library examples postprocessing scripts listed in the file.
e1;e2;...;en        | Run all modules and procedures examples postprocessing scripts matching the semicolon-separated strings.

> The specified file can be automatically generated by the `getList.py`
> Python script in the corresponding language subfolder of the `example`
> folder in the root directory of the project.

> Specific items in a given file can be excluded by
> adding an exclamation mark at the beginning of the line.

> Setting this option to anything other than `none` will lead to running
the Python postprocessing scripts for the corresponding ParaMonte library
examples that have been previously built and run on the command line.

> Frequently, postprocessing scripts merely generate visualizations.
> The visualizations are essential for the ParaMonte library documentation.

> Setting this option to anything other than `none` will lead to running
> the corresponding ParaMonte library examples postprocessing scripts on the command line.

> The examples postprocessing scripts will be executed in the corresponding library example subfolder in the build folder.

**optional**. The default value for `example_postprocessing_list`
is the same as the value set for the option `exam`.

### `fcf`

Specifies any additional compile flags passed to
the Fortran compiler with which the library is built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --fcf "additional_fortran_compiler_flags"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dfcf="additional_fortran_compiler_flags"
    ```

where `additional_fortran_compiler_flags` is string of Fortran compiler flags
separated by semicolon `;` that are added to library default compiler flags.

> The contents of `fcf` is passed to the compiler as is without any changes.

> Any additional preprocessor macros for the compiler can be specified using this option.

> This is low-level build option, mostly useful for the library developers.

**optional**. The default value for `fcf` is an empty string.

### `flf`

Specifies any additional linker flags passed to
the Fortran linker with which the library object files are linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --flf "additional_fortran_linker_flags"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dflf="additional_fortran_linker_flags"
    ```

where `additional_fortran_linker_flags` is string of Fortran linker flags
separated by semicolon `;` that are added to library default linker flags.

> The contents of `flf` is passed to the linker as is without any changes.

> This is low-level build option, mostly useful for the library developers.

**optional**. The default value for `flf` is an empty string.

### `fpp`

Specifies the preprocessing style of the Fortran source files before building the library.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --fpp "fortran_preprocessing_style"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dfpp="fortran_preprocessing_style"
    ```

where `fortran_preprocessing_style` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | No preprocessed source files will be output to the deploy directory.
`default`           | The preprocessed compiler-specific source files will be output to the `fpp` subdirectory.
`generic`           | The preprocessed compiler-agnostic source files will be output to the `fpp` subdirectory.

> This option is relevant to only the Intel and GNU compiler suites.

> When `default` is specified as the value of `fpp`, preprocessed source files may contain
non-portable compiler-specific extensions to the standard Fortran commands and syntax.

> Specifying `default` while building the library with GNU compilers will define macro `GNU_ENABLED=1`.

> Specifying `default` while building the library with Intel compilers will define macro `INTEL_ENABLED=1`.

**optional**. The default value for `fortran_preprocessing_style` depends on the value of `lang` flag.
1.  If the `lang` configuration flag is set to `fortran`, then `fortran_preprocessing_style` is set to,
    1.  `default` if the compiler choice belongs to the family of Intel or GNU compilers.
    2.  `generic` if the compiler choice belongs to other vendors.
2.  If the `lang` configuration flag is set to any possibility other than `fortran`,
    then `fortran_preprocessing_style` is set to `none`.

### `fresh`

Specifies the subdirectories of the build and deployment directories
that must be deleted before starting the new CMake configuration and build.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --fresh "subdirectories_to_delete"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dfresh="subdirectories_to_delete"
    ```

where `subdirectories_to_delete` can be any of the following values.

Value               | Usage
--------------------|------
`all`               | All build and deployment directories for the current library build configuration will be deleted before the new build.
`bdir`              | All build and deployment directories for the current library build configuration will be deleted before the new build.
xx/yy               | The specified `yy` sub-subdirectory in the `xx` subdirectory of the build directory will be deleted before the new build.
filename            | The specified filename in the current build directory will be deleted before the new build.

> This option is simplifies the task of cleaning CMake build directory when library for a given build is being built repeatedly.

> **NOTE**
> Beware that the specified values for all flags (including `--fresh` are sticky and have to be unset for the next builds, if desired.
> For example, the value is set to `all`, any subsequent CMake reconfigurations will use this value even if the flag `--fresh` is not specified
> in the subsequent configurations. This sticky behavior can lead to complete rebuilds of the library, which may be time consuming.
> To change the sticky behavior, one has to reset the value of `--fresh` explicitly.

**optional**. The default value for `subdirectories_to_delete` is `none`.

### `G`

Specifies the CMake build generator.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    -G "cmake_build_generator"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -G "cmake_build_generator"
    ```

where `cmake_build_generator` can be any of the possible CMake build generators.

Value               | Usage
--------------------|------
`NMake Makefiles`   | Generates NMake makefiles.
`MinGW Makefiles`   | Generates makefiles for use with mingw32-make under a Windows command prompt.
`MSYS Makefiles`    | Generates makefiles for use with MSYS (Minimal SYStem) make under the MSYS shell.
`Unix Makefiles`    | Generates standard UNIX makefiles.
`Ninja`             | Generates a `build.ninja` file into the build tree.
other               | Any build generator supported by CMake.

> This is a low-level build setting that is automated by the ParaMonte install scripts.
> This optional flag can be used to enforce a particular CMake build generator when
> the install scripts fail to identify the right build generator for CMake.

> This build flag **must be explicitly set** when CMake executable is directly invoked.

**optional**. Only when specified for the install scripts.
The default value for `cmake_build_generator` is
1.  `Unix Makefiles` on Unix systems (macOS, Linux).
2.  `MSYS Makefiles` on Windows systems within MSYS environments (e.g., MSYS terminal).
4.  `MinGW Makefiles` on Windows systems within Bash MinGW environments (e.g., Git Bash terminal) or Windows CMD terminal.
3.  `NMake Makefiles` on Windows systems within CMD terminal environments.

### `j`

Specifies the number parallel threads to be used by the Make software for parallel library build.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    -j "num_parallel_threads"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -j "num_parallel_threads"
    ```
    or
    ```cmake
    -j
    ```

where `num_parallel_threads` must be a positive integer.

> This option is only relevant to the GNU Make software.

> This option does not impact builds made via Microsoft NMake software as of Jan 2024.

> Specifying a null value when invoking CMake executable binary directly is equivalent to using the maximum available number of threads.

**optional**. The default is the maximum number of parallel threads available on the system.

### `libname`

Specifies the binary name of the output library built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --libname "desired_paramonte_library_name"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dlibname="desired_paramonte_library_name"
    ```

where `desired_paramonte_library_name` points to the
specific name to be used for generating the output library.

**optional**. The default value for `desired_paramonte_library_name` is `libparamonte`.

### `mod`

Specifies the list of desired ParaMonte library modules (all source files beginning with `pm_`)
to compile and add to the final output library file.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --mod "module_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dmod="module_list"
    ```

where `module_list` can be any of the following values.

Value               | Usage
--------------------|------
`all`               | Add and build all modules and their dependencies.
`reset`             | Rebuild the dependencies list for the last specified value of `mod`.
m1;m2;...;mn        | Add and build all modules matching any of the semicolon-separated strings `m1`, ... and their dependencies will be built.

> This option is particularly useful for reducing the library size and
build time by only specifying components are important to the end user.

> To expedite the build process, the generated source file list
> based on the specified value for `mod` is cached within CMake,
> such that subsequent calls to cmake do not regenerate the source
> file list unless the value of `mod` is explicitly changed or set
> to `reset`.

**optional**. The default value for `module_list` is `all` when option `lang` is set to `fortran`
and `pm_sampling` when `lang` is set to any other possible value (all other programming languages).

### `pdt`

Specifies whether the Parameterized Derived Types (PDT) interfaces
of the ParaMonte library should be considered or dropped in the build.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --pdt "pdt_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dpdt="pdt_list"
    ```

where `pdt_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | Deactivates the PDT interfaces within the library  by defining macro `PDT_ENABLED=0`.
`all`               | Activates the PDT interfaces within the library by defining macro `PDT_ENABLED=1`.

> This option primarily exists because of the sparse compiler implementation of
> Parameterized Derived Types (PDT) feature of the Fortran programming language.

> While the PDT interfaces work just fine with the Intel compilers, the current
> implementation status of PDTs in GNU compilers is highly inadequate and leads
> to frequent runtime crashes. This option safely deactivates the PDT interfaces
> and their implementations within the ParaMonte library.

> This is a low level option and mostly useful to the developers of the ParaMonte library.

**optional**. The default is `none`.

### `perfprof`

Specifies whether the library should be compiled for performance profiling.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --perfprof "perfprof_list"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dperfprof="perfprof_list"
    ```

where `perfprof_list` can be any of the following values.

Value               | Usage
--------------------|------
`none`              | Deactivates compiler flags for performance profiling of the library.
`all`               | Activates compiler flags for performance profiling of the library.

> This is a low level option and mostly useful to the developers of the ParaMonte library.

> Performance profiling is currently available only for the GNU Fortran compiler.
> CMake will stop with an error if Intel compilers are specified for performance profiling.
> This option has no effects on the build for all other compiler choices.

**optional**. The default is `none`.

### `purity`

Specifies whether (designated) procedures should be compiled with `pure` or `impure` attribute.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --purity "purity_type"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dpurity="purity_type"
    ```

where `purity_type` can be any of the following values.

Value               | Usage
--------------------|------
`pure`              | Sets the macro `PURE=pure` to add `pure` attribute to the designated procedures.
`impure`            | Sets the macro `PURE=impure` to add `impure` attribute to the designated procedures.

> This is a low level option and mostly useful to the developers of the ParaMonte library.

> Specifying `pure` can potentially lead to faster runtimes by
hinting to the compiler which procedures are `pure` and can be inlined.

> This option must be unset or set to `impure` when `checking` is set to `checked`.

**optional**. The default value for `purity_type` is `pure` if `checking` is set to `nocheck`,
otherwise `impure` if `checking` is set to `checked`.

## TIER-4 ParaMonte library build configuration flags

The ParaMonte TIER-4 build flags set the additional optional arguments that
change the behavior of the library at the lowest level by setting the types.
These build flags are mostly relevant to the ParaMonte library developers or
advanced users who wish to reduce the final library size or supported types.

The values of these build options are automatically set by the CMake scripts.
Changing the default behavior requires careful attention to the consequences.

The Fortran language 2023 supports five intrinsic types:

+   character
+   integer
+   logical
+   complex
+   real

each of which can have a range of possible kind type
parameter values that are specific to the compiler and platform.

The ParaMonte library currently allows building the library for
up to 5 different kind type parameters as specified by the
following variables in the `iso_fortran_env` intrinsic
module of the standard:

Value                   | Usage
------------------------|------
`character_kinds`       | The constant vector of supported character kinds.
`integer_kinds`         | The constant vector of supported integer kinds.
`logical_kinds`         | The constant vector of supported logical kinds.
`real_kinds`            | The constant vector of supported complex/real kinds.

By default, the library is built for all supported kinds that are **fully** supported by the compiler.
This can lead to code bloat and long compilation times if all kinds are activated.
If only specific kinds are desired, the following flags can be used to build the
library only for the selected kinds (specified by their element indices in the
corresponding constant vectors from the `iso_fortran_env` intrinsic module).

> The unicode support within gfortran is an example of a type kind parameter
> that is supported by the compiler but not fully implemented and not functional.
> As such, the type kind parameter corresponding to uniccode is dropped in the library build.

> **WARNING**
> If you specify any of the optional kind-type-parameter-selection command arguments,
> always ensure to minimally specify the default supported kind type parameters by the processor.
> For example, the Fortran standard requires support for at least two `real` kind type parameters.
> Additionally, the standard requires minimal support for at least two `integer` kind type parameters.
> Additionally, the standard requires minimal support for at least one `logical` kind type parameters.

### `ski`

Specifies a list of `character` kind indices
for which the ParaMonte library must be built.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --ski "character_kinds_indices"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dski="character_kinds_indices"
    ```

where `character_kinds_indices` is a semicolon-separated list of
indices of the `character_kinds(:)` constant vector of `iso_fortran_env`
intrinsic module of Fortran.

For example, setting `character_kinds_indices` to `1;3;6` will
build the library for the kinds in `character_kinds([1, 3, 6])`.

> No more than fives indices can be simultaneously specified in every single build, because
> the ParaMonte library is currently configured for builds with up to fives kinds for each type.

> The minimally-required `character` type kind parameter indices for the Intel and GNU compilers are `1`.

**optional**. The default value for `character_kinds_indices`
is all kind numbers fully supported by the compiler.

### `cki`

Specifies a list of `complex` kind indices
for which the ParaMonte library must be built.

Specifies any additional linker flags passed to
the Fortran linker with which the library object files are linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --cki "complex_kinds_indices"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dcki="complex_kinds_indices"
    ```

where `complex_kinds_indices` is a semicolon-separated list of
indices of the `complex_kinds(:)` constant vector of `iso_fortran_env`
intrinsic module of Fortran.

For example, setting `complex_kinds_indices` to `1;3;6` will
build the library for the kinds in `complex_kinds([1, 3, 6])`.

> No more than fives indices can be simultaneously specified in every single build, because
> the ParaMonte library is currently configured for builds with up to 5 kinds for each type.

> The complex kinds can be normally limited to only the two default compiler
> choices for the complex kind without any harm to the library build and run.
> This can lead to noticeable savings in library size and build time.

> The minimally-required `complex` type kind parameter indices for the Intel and GNU compilers are `1;2`.

> **WARNING** Beware some ParaMonte algorithms require compatible `complex` and `real` kind,
> necessitating the use of the same kind type parameters for both data types.

**optional**. The default value for `complex_kinds_indices`
is the kinds specified for the `real` type (below).

### `iki`

Specifies a list of `integer` kind indices
for which the ParaMonte library must be built.

Specifies any additional linker flags passed to
the Fortran linker with which the library object files are linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --iki "integer_kinds_indices"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Diki="integer_kinds_indices"
    ```

where `integer_kinds_indices` is a semicolon-separated list of
indices of the `integer_kinds(:)` constant vector of `iso_fortran_env`
intrinsic module of Fortran.

For example, setting `integer_kinds_indices` to `1;3;6` will
build the library for the kinds in `integer_kinds([1, 3, 6])`.

> No more than fives indices can be simultaneously specified in every single build, because
> the ParaMonte library is currently configured for builds with up to fives kinds for each type.

> If ever specified, always ensure the two default integer kinds supported by
> the Fortran programming language are included in the specified list of indices.

> The minimally-required `integer` type kind parameter indices for the Intel and GNU compilers are `3;4`.

**optional**. The default value for `integer_kinds_indices` is all kinds
supported by the compiler for the C/C++/Fortran programming languages
and only the minimally-required kind type parameters for
all other programming languages.

### `lki`

Specifies a list of `logical` kind indices
for which the ParaMonte library must be built.

Specifies any additional linker flags passed to
the Fortran linker with which the library object files are linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --lki "logical_kinds_indices"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Dlki="logical_kinds_indices"
    ```

where `logical_kinds_indices` is a semicolon-separated list of
indices of the `logical_kinds(:)` constant vector of `iso_fortran_env`
intrinsic module of Fortran.

For example, setting `logical_kinds_indices` to `1;3;6` will
build the library for the kinds in `logical_kinds([1, 3, 6])`.

> The logical kinds can be normally limited to only the default compiler
> choice for the logical kind without any harm to the library build and run.
> This can lead to noticeable savings in library size and build time.

> The minimally-required `logical` type kind parameter indices for the Intel and GNU compilers are `3`.

**optional**. The default value for `logical_kinds_indices` is all kinds
supported by the compiler for the Fortran programming language and
only the minimally-required kind type parameter for
all other programming languages.

### `rki`

Specifies a list of `real` kind indices
for which the ParaMonte library must be built.

Specifies any additional linker flags passed to
the Fortran linker with which the library object files are linked.

+   Usage (with `install.bat` or `install.sh`)
    ```bash
    --rki "real_kinds_indices"
    ```
+   Usage (with `cmake` binary executable)
    ```cmake
    -Drki="real_kinds_indices"
    ```

where `real_kinds_indices` is a semicolon-separated list of
indices of the `real_kinds(:)` constant vector of `iso_fortran_env`
intrinsic module of Fortran.

For example, setting `real_kinds_indices` to `1;3;6` will
build the library for the kinds in `real_kinds([1, 3, 6])`.

> No more than fives indices can be simultaneously specified in every single build, because
> the ParaMonte library is currently configured for builds with up to 5 kinds for each type.

> The real kinds can be normally limited to only the two default compiler
> choices for the real kind without any harm to the library build and run.
> This can lead to noticeable savings in library size and build time.

> The minimally-required `real` type kind parameter indices for the Intel and GNU compilers are `1;2`.

> **WARNING** Beware some ParaMonte algorithms require compatible `complex` and `real` kind,
> necessitating the use of the same kind type parameters for both data types.

**optional**. The default value for `real_kinds_indices` is all kinds
supported by the compiler for the C/C++/Fortran programming languages
and only the minimally-required kind type parameters for
all other programming languages.
