## The ParaMonte library quick start

1.  Clone the ParaMonte repository on your system.
    ```bash
    git clone https://github.com/cdslaborg/paramonte.git
    ```
    Alternatively, you can **optionally** download the **OpenBLAS** and **other ParaMonte submodules** via,
    ```bash
    git clone --recurse-submodules https://github.com/cdslaborg/paramonte.git
    ```

2.  Navigate to the root directory of the repository 
    either in a Windows `CMD` command prompt or Unix `Bash` terminal.
3.  Ensure you have a recent CMake software (`>3.21`), and a recent 
    Intel (`>2021`) or GNU (`>10`) C/Fortran compilers already installed in your terminal.
4.  Using a Windows `CMD` command prompt, type,
    +   For the ParaMonte C library build,
        ```batch
        install.bat --lang c
        ```
    +   For the ParaMonte C++ library build,
        ```batch
        install.bat --lang cpp
        ```
    +   For the ParaMonte Fortran library build,
        ```batch
        install.bat --lang fortran
        ```
    +   For the ParaMonte MATLAB library build,
        ```batch
        install.bat --lang matlab
        ```

    See more on the relevant installation instructions in [install.bat.md](./install.bat.md).    

5.  Using a Unix `Bash` terminal, type,
    +   For the ParaMonte C library build,
        ```bash
        ./install.sh --lang c
        ```
    +   For the ParaMonte C++ library build,
        ```bash
        ./install.sh --lang cpp
        ```
    +   For the ParaMonte Fortran library build,
        ```bash
        ./install.sh --lang fortran
        ```
    +   For the ParaMonte MATLAB library build,
        ```bash
        ./install.sh --lang matlab
        ```
    See more on the relevant installation instructions in [install.sh.md](./install.sh.md).  

4.  For MPI-parallel library builds, ensure you have the relevant MPI library installed on your system.  
    Then, merely add the extra build configuration flag `--par mpi` to the relevant install commands above.  

See 
[install.md](./install.md), 
[install.config.md](./install.config.md), 
[install.bat.md](./install.bat.md), 
[install.sh.md](./install.sh.md), and 
[CMakeLists.md](./CMakeLists.md) 
for more details on library build options and guidelines.  

## Example install commands  

1.  To build for the C, C++, and Fortran programming languages and build and run all their corresponding examples in serial, try, 
    +   On **Windows**,
        ```batch
        install.bat --lang "c;cpp;fortran" --exam all
        ```
    +   On **Linux** and **macOS**,
        ```batch
        ./install.sh --lang "c;cpp;fortran" --exam all
        ```

1.  As of 2024, the following commands can be used to generate 
    minimal `serial`, `openmp` and `mpi` ParaMonte MATLAB binaries.  
    +   On **Windows**,
        ```batch
        install.bat --lang matlab --par "serial;openmp;mpi" --ski 1 --iki "3;4" --lki 3 --rki "1;2" --cki "1;2" --matlabroot "C:\Program Files\MATLAB\R2023a"
        ```
    +   On **Linux**,
        ```batch
        ./install.sh --lang matlab --par serial;openmp;mpi --ski 1 --iki "3;4" --lki 3 --rki "1;2" --cki "1;2" --matlabroot "/usr/local/MATLAB/R2023a"
        ```
    +   On **macOS**,
        ```batch
        ./install.sh --lang matlab --par serial;openmp;mpi --ski 1 --iki "3;4" --lki 3 --rki "1;2" --cki "1;2" --matlabroot "/Applications/MATLAB_R2023b.app"
        ```

See [install.config.md](./install.config.md) for the meaning 
of the flags used and many more possible configuration flags.    