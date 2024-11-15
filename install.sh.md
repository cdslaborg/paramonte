> See [install.md](./install.md) for general installation guidelines.

##  Building via `install.sh` script in Bash terminal (macOS, Linux, MinGW, MSYS2, WSL)

+   The ParMonte library ships with a Bash script file named `install.sh` written in
    the Bash scripting language, located in the root directory of the ParaMonte repository.
+   This script is a convenience wrapper around the lower-level CMake build instructions
    below and automates most installation steps on Unix (Bash) terminals.
+   This is the recommended build mechanism on **macOS** and **Linux** systems.

### Building via `install.sh` script in Bash terminal - quick start

Run the script `install.sh` on a Bash terminal as:

```bash  
./install.sh --lang TARGET_LANAGUAGE
```  

> **Note**
> 
> To redirect the `install.sh` script output to a file named `install.sh.out`, try:
> 
> ```bash  
> ./install.sh --lang TARGET_LANAGUAGE > install.sh.out 2>&1
> ```  
> 
> To redirect and run the `install.sh` script in the background, try:
>
> ```bash  
> ./install.sh --lang TARGET_LANAGUAGE > install.sh.out 2>&1 &
> ```  
>
> To redirect and run the `install.sh` script and disown the process, try:
> 
> ```bash  
> ./install.sh --lang TARGET_LANAGUAGE > install.sh.out 2>&1 & jobs; disown
> ```  
> 
> where you must replace `TARGET_LANAGUAGE` with your choice of programming
> language from which you intend to access the ParaMonte library. 

See the command-line configuration flag descriptions below for a 
list of supported programming languages and other build options.

### Building via `install.sh` script in Bash terminal - prerequisites

You will need the following components installed
on your PC to successfully run `install.sh` script.

1.  A recent version (`>2018`) installation of **CMake** build generator software.  
    This application can be freely [downloaded and installed](https://cmake.org/download/).

2.  A recent **C and Fortran compilers** minimally supporting Fortran 2008.  

    + GNU and Intel compilers are popular compiler choices for **Linux** systems.  

    + GNU and NAG compilers are popular compiler choices for **macOS** systems.  

        Although NAG is an excellent Fortran compiler choice, it is untested.

    +   Two popular compiler choices for **WSL**, **MSYS2**, and **MinGW** Windows 
        environments for systems are GNU and Intel compilers.  

    +   Building via **Intel compilers** requires a recent version (`>2023`) 
        installation of Intel OneAPI and HPC Toolkits on your system.  
        The Intel OneAPI can be freely downloaded and installed by anyone.

    +   Building via GNU compilers requires a recent version (`>10`) 
        installation of GNU C and Fortran compilers.  

    +   Building for **MPI-parallelism** requires one of the following MPI libraries.

        +   The Intel MPI library with Intel compilers.
            The Intel MPI library can be downloaded as part of the 
            Intel OneAPI + HPC compiler suite from the Intel website.
        +   The MPICH MPI library with GNU compilers. 
            For download instructions visit the MPICH library website.
            You can also install the MPICH library via your preferred package managers.
            For example, 
            +   on **Linux Ubuntu** operating systems using the APT package manager, 
                you can try the following command in a Bash terminal to install MPICH MPI library:
                ```bash
                sudo apt install mpich
                ```
            +   on **macOS (Darwin)** operating systems using the `brew` package manager, 
                you can try the following command in a ZSH or Bash terminal to install MPICH MPI library:
                ```bash
                brew install mpich && brew link mpich
                ```
        +   The Open-MPI library with GNU compilers.
            For download instructions visit the Open-MPI library website.
            You can also install the Open-MPI library via your preferred package managers.
            For example, 
            +   on **Linux Ubuntu** operating systems using the APT package manager, 
                you can try the following command in a Bash terminal to install Open-MPI MPI library:
                ```bash
                sudo apt install openmpi-bin openmpi-common libopenmpi-dev
                ```
            +   on **macOS (Darwin)** operating systems using the `brew` package manager, 
                you can try the following command in a ZSH or Bash terminal to install Open-MPI MPI library:
                ```bash
                brew install openmpi && brew link openmpi
                ```
        +   The Microsoft MPI library with GNU compilers (on Windows).
            The ParaMonte library has never been tested against Microsoft MPI library.
            Microsoft claims full compatibility of their MPI library with the MPI standard.
            Yet, it is unclear whether their library has any Fortran bindings.
            If you have any luck building ParaMonte and linking it against 
            Microsoft MPI library, please share your success with us at
            [ParaMonte GitHub Disussions Page](https://github.com/cdslaborg/paramonte/discussions).

3.  Once all the above components are installed, open a Windows CMD terminal that
    recognizes all applications installed above (that is, all applications can
    be found in the `PATH` environment variable).

6.  Navigate to the folder containing the ParaMonte repository on your system.

7.  Type `install.sh` with the desired build configuration flags. Example:

    ```batch
    ./install.sh --lang c
    ```

    ```bash
    ./install.sh --lang c --build release --par mpi --checking nocheck
    ```

    ```bash
    ./install.sh --lang c --build release --par "mpi;omp" --checking nocheck
    ```

    Check out [install.config.md](./install.config.md) file for a full list of build configuration flags.
