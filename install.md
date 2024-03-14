## The ParaMonte library build mechanisms

There are three ways to build the ParaMonte library:

1.  Building via the [install.sh](https://github.com/cdslaborg/paramonte/blob/main/install.sh) Bash script located in the root 
    directory of the [ParaMonte GitHub repository](https://github.com/cdslaborg/paramonte) in the in a **Unix Bash terminal**.  
    See the relevant installation instructions in [install.sh.md](./install.sh.md).  
    This approach works on all platforms that support Bash terminals including:  

    +   **Git Bash** or other **MinGW** terminals on **Windows**
    +   **MSYS2** on **Windows**
    +   **WSL** on **Windows**
    +   **macOS**
    +   **Linux** 

    See the relevant installation instructions in [install.sh.md](./install.bat.md).  

2.  Building via the [install.bat](https://github.com/cdslaborg/paramonte/blob/main/install.bat) Batch script located in the root 
    directory of the [ParaMonte GitHub repository](https://github.com/cdslaborg/paramonte) in the in a **Windows CMD terminal**.  
    See the relevant installation instructions in [install.bat.md](./install.bat.md).    

3.  Building via the CMake script directly in any Windows or Unix terminal that supports CMake.  
    See the relevant installation instructions in [CMakeLists.md](./CMakeLists.md).  

The Bash and Batch install scripts in the first two build mechanisms above are merely 
convenient wrappers around the lower-level CMake scripts in the third building mechanism.  
