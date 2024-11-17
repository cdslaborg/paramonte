## Using the CMake build-generator software in any Linux, macOS, or Windows terminal.

The ParMonte library ships with an extensive set of CMake scripts that build the library. 
While it is possible to build the library via CMake, the build process involves several
steps that must be carried out manually by the user, in addition to what is normally 
required for builds via the [install.bat](./install.bat.md) 
and [install.sh](./install.sh.md) installation scripts.  

For more information on automated build mechanisms, see the contents of [install.md](./install.md).

### The ParaMonte library build instructions via CMake

#### Windows CMD shell

On Windows platforms using Visual Studio NMake within a Windows CMD 
terminal (with Intel version `>2021` or GNU Fortran version `>10` compilers),

+   Ensure you have a recent (`>19`) Microsoft Visual Studio
    installed on your system (the free Community Edition is enough).

+   Open a Windows command prompt that recognizes the Microsoft NMake application
    along with the compiler you intend to use (e.g., Intel or GNU compilers).

    +   If you have the Intel compilers installed and you intend to use it,
        we strongly recommend you open the Intel OneAPI Command
        Prompt that automatically ships with Intel products.
        The Intel CMD is automatically preloaded with
        compiler/library definitions and paths.

    +   If you intend to configure with GNU or other compilers, you can open
        any CMD terminal that minimally recognizes the Microsoft NMake application
        and your compiler of interest. We recommend using the CMD
        terminal that ships automatically with Microsoft Visual Studio.

+   Navigate to the root directory of the ParaMonte library
    (where the library LICENSE file exists).

+   Create a `build` directory within the project root
    and change the current directory to it via,
    ```batch
    mkdir build && cd build
    ```

+   Call the CMake executable via the following command,
    ```batch
    cmake .. -G "NMake Makefiles" <options>
    ```
    where `<options>` must be removed or replaced with any
    of the configuration options described later below.

+   Once successfully configured, build the library by calling Microsoft `nmake`,
    ```batch
    nmake
    ```

+   Install the built library in the `lib` subdirectory within the build directory via,
    ```batch
    nmake install
    ```

In sum, you can use the following all-in-one minimal command to build the library,
```batch
mkdir build & cd build && cmake -G "NMake Makefiles" .. && nmake && nmake install
```


Once the library is built, you can also further take additional optional steps:

+   Deploy the library to the prespecified deployment folder via,
    ```batch
    nmake deploy
    ```

+   Test the library (if enabled at CMake configuration time) via,
    ```batch
    nmake test
    ```

+   Build and run the library examples (if enabled at CMake configuration time) via,
    ```batch
    nmake example
    ```

+   Build and run the library benchmarks (if enabled at CMake configuration time) via,
    ```batch
    nmake benchmark
    ```

#### Unix Bash shell

On **Windows**, **WSL**, or **Unix** (**Linux**/**Darwin**/..) platforms using a POSIX-style terminal
(e.g., **Bash**, **MSYS2**, **Windows Git Bash**, ...) along with GNU (**MinGW**) Make application,

+   Open a POSIX-compatible terminal that minimally recognizes,

    +   A recent version of CMake (`>3.16`),
    +   The GNU Make application (Unix Makefiles or MinGW Makefiles),
    +   A Fortran 2008-compliant compiler (e.g., gfortran `>10`, ifort `>2021`, ifx `>2024`)

+   Navigate to the root directory of the ParaMonte library
    (where the library LICENSE file exists).

+   Create a `build` directory within the project root
    and change the current directory to it via,
    ```bash
    mkdir build && cd build
    ```

+   Call CMake to build the application,

    +   On Windows MinGW terminals (e.g., Git Bash), try,
    ```bash
    cmake -G "MinGW Makefiles" .. -Dfc="$(command -v gfortran.exe)"
    ```

    +   On Windows MSYS2 terminals, try,
        ```bash
        cmake -G "MSYS2 Makefiles" .. -Dfc="$(command -v gfortran.exe)"
        ```

    +   On all other POSIX-compatible terminals (including WSL), try,
        ```bash
        cmake -G "Unix Makefiles" .. -Dfc="$(command -v gfortran)"
        ```

    where you can change `gfortran` or `gfortran.exe` with any other
    Fortran compiler that you intend to use. Feel free to supply any
    other build configuration flags described later below.

+   Once successfully configured, build the library by calling the GNU `make`,
    ```bash
    make
    ```

    +   If using a Windows MSYS2 terminal, install the GNU `make` software.
    +   If using a Windows MinGW terminal, install the GNU `mingw32-make` software.
        If the `make` software name is `mingw32-make` within the MinGW terminal, either use this name
        or change the binary name `mingw32-make` to `make` before calling `make`.

+   Install the built library in the `lib` subdirectory within the build directory via,
    ```bash
    make install
    ```

In sum, you can use the following all-in-one command within MinGW environments,
```bash
mkdir build; cd build && cmake -G "MinGW Makefiles" .. -Dfc="$(command -v gfortran.exe)" && make && make install
```


In sum, you can use the following all-in-one command within MSYS2 environments,
```bash
mkdir build; cd build && cmake -G "MSYS2 Makefiles" .. -Dfc="$(command -v gfortran.exe)" && make && make install
```

In sum, you can use the following sequence of commands within all other Unix-compatible environments,
```bash
mkdir build; cd build && cmake -G "Unix Makefiles" .. -Dfc="$(command -v gfortran)" && make && make install
```


Once the library is built, you can also further take additional optional steps:

+   Deploy the library to the prespecified deployment folder via,
    ```bash
    make deploy
    ```

+   Test the library (if enabled at CMake configuration time) via,
    ```bash
    make test
    ```

+   Build and run the library examples (if enabled at CMake configuration time) via,
    ```bash
    make example
    ```

+   Build and run the library benchmarks (if enabled at CMake configuration time) via,
    ```bash
    make benchmark
    ```
