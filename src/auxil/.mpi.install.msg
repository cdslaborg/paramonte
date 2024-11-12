The MPI runtime libraries are required ONLY for the parallel ParaMonte simulations.
The MPI-parallel version of the ParaMonte library requires
one of the following MPI runtime libraries:

    1)  Intel MPI runtime library (recommended and available only for Windows and Linux)

        A)  Download the Intel MPI runtime library from the Intel website and install it.
        B)  Ensure the mpiexec MPI executable binary is available on the command line in the terminal.
            I)  On Windows systems, you may have difficulty finding mpiexec in Windows terminals.
                For this reason, Intel provides its dedicated Intel OneAPI command-prompt
                on Windows which has all environmental paths predefined within the shell.
                This Batch terminal may not ship with MPI runtime library installations.
                If so, you can directly run the Intel Windows `vars.bat` Batch script within a
                compatible Windows shell to define all environmental variables for MPI-parallel simulations.
                By default, this script is frequently (but not always) stored in one of the following path:

                    "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
                    "C:\Program Files (x86)\Intel\oneAPI\mpi\latest\env\vars.bat"

            II) On Linux systems, you can also directly run the Intel `vars.sh` Bash script within a
                compatible terminal to define all necessary environmental variables for MPI-parallel simulations.
                By default, this script is frequently (but not always) stored in the following installation path:

                    "/opt/intel/oneapi/mpi/latest/env/vars.sh"

    2)  MPICH MPI runtime library (available only for Linux and macOS/Darwin)

        This can be downloaded and installed via a package manager.
        On Linux Ubuntu systems, you can use the shell installation command:

            sudo apt install mpich

        On macOS/Darwin systems, you can use the Homebrew application to install MPICH:

            brew install mpich && brew link mpich

            If you do not have Homebrew installed on your system, you can do so by
            typing the following command in your (Bash/zsh) shell terminal:

                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    3)  OpenMPI runtime library (available only for Linux and macOS/Darwin)

        This can be downloaded and installed via a package manager.
        On Linux Ubuntu systems, you can use the shell installation command:

            sudo apt install openmpi-bin openmpi-common libopenmpi-dev

        On macOS/Darwin systems, you can use the Homebrew application to install OpenMPI:

            brew install openmpi && brew link openmpi

            If you do not have Homebrew installed on your system, you can do so by
            typing the following command in your (Bash/zsh) shell terminal:

                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

For more information on running MPI-parallel simulations, visit:

    <a href = "www.cdslab.org/pm">www.cdslab.org/pm</a>