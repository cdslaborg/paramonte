%>  \brief
%>  Return ``true`` if and only if the application has been invoked via
%>  an Intel MPI ``mpiexec`` launcher, otherwise, return ``false``.<br>
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
%>  </ol>
%>
%>  Therefore, this function attempts to detect the above runtime environment variables.<br>
%>
%>  \return
%>  ``itis``    :   The output scalar MATLAB ``logical`` that is ``true``
%>                  if and only if the application has been invoked via an
%>                  Intel MPI ``mpiexec`` launcher, otherwise, return ``false``.<br>
%>
%>  \interface{isimpi}
%>  \code{.m}
%>
%>      itis = pm.lib.mpi.runtime.isimpi()
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
%>  \example{isimpi}
%>  \include{lineno} example/lib/mpi/runtime/isimpi/main.m
%>  \output{isimpi}
%>  \include{lineno} example/lib/mpi/runtime/isimpi/main.out.m
%>
%>  \final{isimpi}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, 12:10 AM Wednesday, November 13, 2024, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function itis = isimpi()

    brand = getenv("I_MPI_INFO_BRAND");
    state = getenv("I_MPI_INFO_STATE");
    itis = ~isempty(brand) || ~isempty(state);

end