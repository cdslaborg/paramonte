%>  \brief
%>  Return the MPI image count and the current image ID (e.g., MPI rank + 1)
%>  if the application has been launched via the MPICH MPI ``mpiexec`` launcher.<br>
%>
%>  \details
%>  This function attempts to detect the invocation of the current MATLAB
%>  session via an MPICH MPI ``mpiexec`` binary MPI launcher instance.<br>
%>  The runtime detection is achieved by checking for runtime environment
%>  variables that the MPICH MPI ``mpiexec`` binary executable defines upon launch.<br>
%>
%>  Specifically, MPICH MPI provides a long informal list of
%>  environment variables that are defined on every MPI process.<br>
%>  These variables may not be defined on all platforms.<br>
%>  Some are only defined at runtime, while others are defined at compile-time.<br>
%>
%>  <ol>
%>      <li>    On **Unix systems**, the best approach to investigate the MPICH-specific runtime variables
%>              is to launch the MPICH ``mpiexec`` binary executable with the ``printenv`` Unix command.<br>
%>              ```bash
%>                  mpiexec -np 1 printenv
%>              ```
%>              To see the compile-time variables and compare them with the runtime list from the above, try:<br>
%>              ```bash
%>                  printenv
%>              ```
%>      <li>    On **Windows systems**, the best approach to investigate the MPICH-specific runtime variables
%>              is to launch the MPICH ``mpiexec`` binary executable with the ``set`` Windows batch command.<br>
%>              ```batch
%>                  mpiexec -localonly -n 1 cmd /k set
%>              ```
%>              To see the compile-time variables and compare them with the runtime list from the above, try:<br>
%>              ```batch
%>                  set
%>              ```
%>  </ol>
%>
%>  The most persistent environment variables for MPICH ``mpiexec`` binary appear to be the following:<br>
%>
%>  <ol>
%>      <li>    ``PMI_RANK``            :   The rank of the current process in the current MPI communication world.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and MPICH MPI libraries.<br>
%>      <li>    ``PMI_SIZE``            :   The size (i.e., number of images/processes) of the current MPI communication world.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and MPICH MPI libraries.<br>
%>      <li>    ``MPI_LOCALRANKID``     :   The **local** rank of the current process relative to the local host.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and MPICH MPI libraries.<br>
%>      <li>    ``MPI_LOCALNRANKS``     :   The size (i.e., number of images/processes) of the current MPI communication on the local host.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and MPICH MPI libraries.<br>
%>  </ol>
%>
%>  Therefore, this function first attempts to detect the Hydra process manager.<br>
%>  If identified, the function will then attempt to identify the MPICH ``mpiexec`` runtime variables.<br>
%>
%>  \return
%>  ``nproc``   :   The output scalar MATLAB non-negative whole number containing
%>                  the number of MPI images in the current invocation of the
%>                  MPICH MPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible MPICH MPI library or failed to infer the image count.<br>
%>  ``rankp1``  :   The output scalar MATLAB non-negative whole number containing
%>                  the image ID (e.g., MPI rank + 1) of the current MPI image in
%>                  the current invocation of the MPICH MPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible MPI library or failed to infer the image count.<br>
%>                  Note that the image ID always starts at ``1``, unlike the MPI rank.<br>
%>                  The argument ``rankp1`` stands for ``rank + 1``.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible MPICH MPI library or failed to infer the image count.<br>
%>
%>  \interface{mmpi}
%>  \code{.m}
%>
%>      [nproc, rankp1] = pm.lib.mpi.runtime.mmpi()
%>
%>  \endcode
%>
%>  \see
%>  [pm.lib.mpi.runtime.mmpi()](@ref mmpi)<br>
%>  [pm.lib.mpi.runtime.impi()](@ref impi)<br>
%>  [pm.lib.mpi.runtime.ompi()](@ref ompi)<br>
%>  [pm.lib.mpi.runtime.hydra()](@ref hydra)<br>
%>  [pm.lib.mpi.runtime.nproc()](@ref nproc)<br>
%>  [pm.lib.mpi.runtime.rankp1()](@ref rankp1)<br>
%>  [pm.lib.mpi.runtime.isimpi()](@ref isimpi)<br>
%>  [pm.lib.mpi.runtime.detect()](@ref detect)<br>
%>
%>  \example{mmpi}
%>  \include{lineno} example/lib/mpi/runtime/mmpi/main.m
%>  \output{mmpi}
%>  \include{lineno} example/lib/mpi/runtime/mmpi/main.out.m
%>
%>  \final{mmpi}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, 12:10 AM Wednesday, November 13, 2024, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function [nproc, rankp1] = mmpi()

    %%%%
    %%%% Ensure the consistency of environment variables and non-Intel identity of ``mpiexec``.
    %%%%

    if  pm.lib.mpi.runtime.isimpi()% || nproc == 0 || rankp1 == 0
        rankp1 = 0;
        nproc = 0;
    else
        [nproc, rankp1] = pm.lib.mpi.runtime.hydra();
    end

end