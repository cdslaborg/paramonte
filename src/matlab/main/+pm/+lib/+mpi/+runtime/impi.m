%>  \brief
%>  Return the MPI image count and the current image ID (e.g., MPI rank + 1)
%>  if the application has been launched via the Intel MPI ``mpiexec`` launcher.<br>
%>
%>  \details
%>  This function attempts to detect the invocation of the current MATLAB
%>  session via an Intel MPI ``mpiexec`` binary MPI launcher instance.<br>
%>  The runtime detection is achieved by checking for runtime environment
%>  variables that the Intel MPI ``mpiexec`` binary executable defines upon launch.<br>
%>
%>  Specifically, Intel MPI provides a long informal list of
%>  environment variables that are defined on every MPI process.<br>
%>  These variables may not be defined on all platforms.<br>
%>  Some are only defined at runtime, while others are defined at compile-time.<br>
%>
%>  <ol>
%>      <li>    On **Unix systems**, the best approach to investigate the Intel-specific runtime variables
%>              is to launch the Intel ``mpiexec`` binary executable with the ``printenv`` Unix command.<br>
%>              \code{.sh}
%>                  mpiexec -np 1 printenv
%>              \endcode
%>              To see the compile-time variables and compare them with the runtime list from the above, try:<br>
%>              \code{.sh}
%>                  printenv
%>              \endcode
%>      <li>    On **Windows systems**, the best approach to investigate the Intel-specific runtime variables
%>              is to launch the Intel ``mpiexec`` binary executable with the ``set`` Windows batch command.<br>
%>              \code{.bat}
%>                  mpiexec -localonly -n 1 cmd /k set
%>              \endcode
%>              To see the compile-time variables and compare them with the runtime list from the above, try:<br>
%>              \code{.bat}
%>                  set
%>              \endcode
%>  </ol>
%>
%>  Note that all Intel-MPI-specific environment variables are prefixed with ``I_MPI_``.<br>
%>  The most persistent environment variables for Intel ``mpiexec`` binary appear to be the following:<br>
%>
%>  <ol>
%>      <li>    ``I_MPI_INFO_STATE``    :   The status code for the Intel MPI runtime (e.g., ``0``).<br>
%>      <li>    ``I_MPI_INFO_BRAND``    :   The processor brand (e.g., ``12th Gen Intel(R) Core(TM) i9-12950HX``).<br>
%>      <li>    ``PMI_RANK``            :   The rank of the current process in the current MPI communication world.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and Intel MPI libraries.<br>
%>      <li>    ``PMI_SIZE``            :   The size (i.e., number of images/processes) of the current MPI communication world.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and Intel MPI libraries.<br>
%>      <li>    ``MPI_LOCALRANKID``     :   The **local** rank of the current process relative to the local host.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and Intel MPI libraries.<br>
%>      <li>    ``MPI_LOCALNRANKS``     :   The size (i.e., number of images/processes) of the current MPI communication on the local host.<br>
%>                                          This environment variable is common to all Hydra-based MPI libraries, such as MPICH and Intel MPI libraries.<br>
%>  </ol>
%>
%>  Therefore, this function first attempts to detect the Hydra process manager.<br>
%>  If identified, the function will then attempt to identify the Intel ``mpiexec`` runtime variables.<br>
%>
%>  \return
%>  ``nproc``   :   The output scalar MATLAB non-negative whole number containing
%>                  the number of MPI images in the current invocation of the
%>                  Intel MPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible Intel MPI library or failed to infer the image count.<br>
%>  ``rankp1``  :   The output scalar MATLAB non-negative whole number containing
%>                  the image ID (e.g., MPI rank + 1) of the current MPI image in
%>                  the current invocation of the Intel MPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible MPI library or failed to infer the image count.<br>
%>                  Note that the image ID always starts at ``1``, unlike the MPI rank.<br>
%>                  The argument ``rankp1`` stands for ``rank + 1``.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible Intel MPI library or failed to infer the image count.<br>
%>
%>  \interface{impi}
%>  \code{.m}
%>
%>      [nproc, rankp1] = pm.lib.mpi.runtime.impi()
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
%>  \example{impi}
%>  \include{lineno} example/lib/mpi/runtime/impi/main.m
%>  \output{impi}
%>  \include{lineno} example/lib/mpi/runtime/impi/main.out.m
%>
%>  \final{impi}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, 12:10 AM Wednesday, November 13, 2024, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function [nproc, rankp1] = impi()

    %%%%
    %%%% Ensure the consistency of environment variables and Intel identity of ``mpiexec``.
    %%%%

    if ~pm.lib.mpi.runtime.isimpi()% || nproc == 0 || rankp1 == 0
        rankp1 = 0;
        nproc = 0;
    else
        [nproc, rankp1] = pm.lib.mpi.runtime.hydra();
    end

end