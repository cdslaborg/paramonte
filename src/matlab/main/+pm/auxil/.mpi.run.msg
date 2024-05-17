On Linux/macOS systems, launch the MPI-parallel ParaMonte
simulation on the shell terminal via the following command:

    mpiexec -n 3 matlab -batch "main"

where,

    1.  ``mpiexec`` is the MPI launcher application.
    2.  ``-n 3`` specified three processes to be used for the MPI-parallel simulation.
    3.  ``matlab -batch`` launches the LANG_NAME software per MPI process on the command line.
    5.  ``main`` is the name of your script containing the ParaMonte library task to run in parallel.
        Replace the script name with whatever the actual name of your script file is.
        Notice the MATLAB filename `main` appears without the file extension.
        Notice the MATLAB filename is enclosed with quotation marks.

On Windows systems using the Intel MPI library, an additional ``-localonly`` flag is required:

    mpiexec -localonly -n 3 matlab -batch "main"

where,

    1.  ``-localonly`` indicates a parallel simulation on only a single node.
        This flag will obviate the need for MPI library credentials registration.
        For more information, visit the Intel MPI library website.