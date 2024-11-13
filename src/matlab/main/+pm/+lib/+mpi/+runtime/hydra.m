%>  \brief
%>  Return the MPI image count and the current image ID (e.g., MPI rank + 1)
%>  if the application has been launched via a Hydra-based MPI ``mpiexec`` launcher.<br>
%>
%>  \details
%>  This function attempts to detect the invocation of the current MATLAB
%>  session via a Hydra-based MPI ``mpiexec`` binary MPI launcher instance.<br>
%>  The runtime detection is achieved by checking for runtime environment
%>  variables that a Hydra-based ``mpiexec`` binary executable defines upon launch.<br>
%>
%>  Specifically, a Hydra-basedMPI provides an informal list of
%>  environment variables that are defined on every MPI process.<br>
%>  These variables may not be defined on all platforms.<br>
%>  Some are only defined at runtime, while others are defined at compile-time.<br>
%>
%>  <ol>
%>      <li>    On **Unix systems**, the best approach to investigate the Hydra-MPI-specific runtime variables
%>              is to launch a Hydra-based ``mpiexec`` binary executable with the ``printenv`` Unix command.<br>
%>              \code{.sh}
%>                  mpiexec -np 1 printenv
%>              \endcode
%>              To see the compile-time variables and compare them with the runtime list from the above, try:<br>
%>              \code{.sh}
%>                  printenv
%>              \endcode
%>      <li>    On **Windows systems**, the best approach to investigate the Hydra-specific runtime variables
%>              is to launch a Hydra-based ``mpiexec`` binary executable with the ``set`` Windows batch command.<br>
%>              \code{.bat}
%>                  mpiexec -n 1 cmd /k set
%>              \endcode
%>              To see the compile-time variables and compare them with the runtime list from the above, try:<br>
%>              \code{.bat}
%>                  set
%>              \endcode
%>              Note that as of 2024, the only ParaMonte-compatible Hydra-based
%>              MPI library on Windows is the Intel MPI library.<br>
%>  </ol>
%>
%>  The most persistent environment variables for Hydra-based ``mpiexec`` binary appear to be the following:<br>
%>
%>  <ol>
%>      <li>    ``PMI_RANK``                    :   The rank of the current process in the current MPI communication world.<br>
%>      <li>    ``PMI_SIZE``                    :   The size (i.e., number of images/processes) of the current MPI communication world.<br>
%>      <li>    ``MPI_LOCALRANKID``             :   The **local** rank of the current process relative to the local host.<br>
%>                                                  This environment variable may not be defined on all platforms (e.g., Windows) by all Hydra-based MPI libraries.<br>
%>      <li>    ``MPI_LOCALNRANKS``             :   The size (i.e., number of images/processes) of the current MPI communication on the local host.<br>
%>                                                  This environment variable may not be defined on all platforms (e.g., Windows) by all Hydra-based MPI libraries.<br>
%>  </ol>
%>
%>  \return
%>  ``nproc``   :   The output scalar MATLAB non-negative whole number containing
%>                  the number of MPI images in the current invocation of a
%>                  Hydra-based MPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible Hydra-based MPI library or failed to infer the image count.<br>
%>  ``rankp1``  :   The output scalar MATLAB non-negative whole number containing
%>                  the image ID (e.g., MPI rank + 1) of the current MPI image in
%>                  the current invocation of a Hydra-based ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible MPI library or failed to infer the image count.<br>
%>                  Note that the image ID always starts at ``1``, unlike the MPI rank.<br>
%>                  The argument ``rankp1`` stands for ``rank + 1``.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible Hydra-based MPI library or failed to infer the image count.<br>
%>
%>  \interface{hydra}
%>  \code{.m}
%>
%>      [nproc, rankp1] = pm.lib.mpi.runtime.hydra()
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
%>  \example{hydra}
%>  \include{lineno} example/lib/mpi/runtime/hydra/main.m
%>  \output{hydra}
%>  \include{lineno} example/lib/mpi/runtime/hydra/main.out.m
%>
%>  \final{hydra}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, 12:10 AM Wednesday, November 13, 2024, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function [nproc, rankp1] = hydra()

    for varname = ["PMI_SIZE", "MPI_LOCALNRANKS"]
        nproc = getenv(varname);
        if ~isempty(nproc)
            nproc = str2double(nproc);
            break;
        else
            nproc = [];
        end
    end
    if  isempty(nproc)
        nproc = 0;
    end

    for varname = ["PMI_RANK", "MPI_LOCALRANKID"]
        rankp1 = getenv(varname);
        if ~isempty(rankp1)
            rankp1 = str2double(rankp1) + 1;
            break;
        else
            rankp1 = [];
        end
    end
    if  isempty(rankp1)
        rankp1 = 0;
    end

end