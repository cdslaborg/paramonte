> See [install.md](./install.md) for general installation guidelines.

> On Windows platforms, the length of the default library build path 
> is near the maximum value (`250` characters) allowed by CMake software.
> This can lead to build failures, particularly for library examples and benchmarks.
> There are two solutions to this limitation on the Windows OS,
> 1.    Place the ParaMonte library repository near the root Drive on the system,
>       for example, in the `C:\`, `D:\` or any other system drive available.
>       Placing the library in the drive root directory will shorten the length 
>       of the full build path and will likely resolve the build failures.
> 2.    Specify a custom short path for the build directory via the 
>       installation configuration flag [--bdir](./install.config.md#bdir).

##  Building via `install.bat` script in Windows CMD terminal

+   The ParMonte library ships with a Windows Batch script file named `install.bat`
    written in Windows CMD Batch scripting language, located in the root directory
    of the ParaMonte repository.
+   This script is a convenience wrapper around the lower-level CMake build
    instructions below and automates most installation steps
    on Windows systems using the Windows CMD terminal.
+   This script is the recommended build mechanism on **Windows** systems.

### Building via `install.bat` script in Windows CMD terminal - quick start

Run the script `install.bat` on an Intel-aware Command Prompt for Windows as:

```batch  
install.bat --lang TARGET_LANAGUAGE
```  

and to redirect the `install.bat` script output to a file named `install.bat.out`, try:

```batch  
install.bat --lang TARGET_LANAGUAGE > install.bat.out 2>&1
```  

where you must replace `TARGET_LANAGUAGE` with your choice of programming
language from which you intend to access the ParaMonte library. See the
command-line configuration flag descriptions below for a list of
supported programming languages and other build options.

### Building via `install.bat` script in Windows CMD terminal - prerequisites

You will need the following components installed
on your PC to successfully run `install.bat` script.

1.  A recent version (`>2018`) installation of **CMake** build generator software.  
    This application can be freely [downloaded and installed](https://cmake.org/download/).

2.  A recent **C + Fortran compilers** minimally supporting Fortran 2008.  
    Two popular compiler choices for Windows are GNU and Intel compilers.  
    NAG is another excellent Fortran compiler choice, although untested.

    +   Using **Intel Compiler Collection** on Windows CMD

        a.  First, install a recent version of Microsoft Visual Studio (VS) (`>2020`).  
            The Community Edition version of Microsoft Visual Studio can be freely 
            [downloaded](https://visualstudio.microsoft.com/vs/community/) and installed.  
            Ensure C++ development tools are selected for installation as they are
            required to integrate Visual Studio with Intel compilers.

        b.  A recent version (`>2023`) installation of Intel OneAPI Base Toolkit.  
            This application can be freely downloaded and installed by anyone.

        c.  A recent version (`>2023`) installation of Intel OneAPI HPC Toolkit.  
            This application can be freely downloaded and installed by anyone.

        d.  Once all the above components are installed, open the
            Intel Command Prompt from the Windows Start Menu.  
            Below is an illustrative image.  
            ![intel-command-prompt-windows.png](https://raw.githubusercontent.com/cdslaborg/paramonte/refs/heads/main/img/intel-command-prompt-windows.png).

    +   Using **GNU Compiler Collection** on Windows CMD

        a.  A recent version (`>10`) installation of GNU C/C++/Fortran compilers.  
            The [quickstart-fortran](https://github.com/LKedward/quickstart-fortran/releases)
            offers an excellent packing of these tools.

        b.  A recent GNU-compatible MPI library, only if MPI parallelism is desired.  
            Other than the Intel MPI library that ships and integrates with Intel compilers,
            the most promising alternative on the Windows Operating System seems to be the
            [Microsoft MPI library](https://github.com/microsoft/Microsoft-MPI).  
            The ParaMonte library build has not been tested with the Microsoft MPI library.  
            Only builds with Intel MPI library and compilers are currently tested.

        c.  Once all the above components are installed, open a Windows CMD terminal that
            recognizes all applications installed above (that is, all applications can
            be found in the `PATH` environment variable).  
            Two excellent choices are the **Intel Command Prompt** that 
            can be opened from the Windows Start Menu, as illustrated below.  
            ![intel-command-prompt-windows.png](https://raw.githubusercontent.com/cdslaborg/paramonte/refs/heads/main/img/intel-command-prompt-windows.png)
            and the **x64 Native Tools Command Prompt for VS** that 
            can be opened from the Windows Start Menu, as illustrated below.  
            ![msvs-command-prompt-windows.png](https://raw.githubusercontent.com/cdslaborg/paramonte/refs/heads/main/img/msvs-command-prompt-windows.png)

3.  Navigate to the folder containing the ParaMonte repository on your system.

4.  Type `install.bat` with the desired build configuration flags. Example:

    ```batch
    install.bat --lang c
    ```

    ```batch
    install.bat --lang c --build release --par mpi --checking nocheck
    ```

    ```batch
    install.bat --lang c --build release --par "mpi;omp" --checking nocheck
    ```

    Check out [install.config.md](./install.config.md) file for a full list of build configuration flags.
