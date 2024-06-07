## The ParaMonte library quick start

1.  Clone the ParaMonte repository on your system.
    ```bash
    git clone https://github.com/cdslaborg/paramonte.git
    ```
    Alternatively, you can **optionally** download the **OpenBLAS** and **other ParaMonte submodules** via,
    ```bash
    git clone --recurse-submodules https://github.com/cdslaborg/paramonte.git
    ```

2.  Navigate to the root directory of the repository either in a Windows `CMD` command prompt or Unix `Bash` terminal.
3.  Ensure you have a recent CMake software (`>3.16`), and a recent Intel (`>2021`) or GNU (`>10`) C/Fortran compilers already installed in your terminal.
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
        install.sh --lang fortran
        ```
    +   For the ParaMonte MATLAB library build,
        ```bash
        install.sh --lang matlab
        ```
    See more on the relevant installation instructions in [install.sh.md](./install.sh.md).    

See [install.md](./install.md), [install.config.md](./install.config.md), [install.bat.md](./install.bat.md), [install.sh.md](./install.sh.md), and [CMakeLists.md](./CMakeLists.md) 
for more library build options and guidelines.
