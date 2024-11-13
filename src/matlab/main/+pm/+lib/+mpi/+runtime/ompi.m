%>  \brief
%>  Return the MPI image count and the current image ID (e.g., MPI rank + 1)
%>  if the application has been launched via the OpenMPI ``mpiexec`` launcher.<br>
%>
%>  \details
%>  This function attempts to detect the invocation of the current MATLAB
%>  session via an OpenMPI ``mpiexec`` binary MPI launcher instance.<br>
%>  The runtime detection is achieved by checking for runtime environment
%>  variables that the OpenMPI ``mpiexec`` binary executable defines upon launch.<br>
%>
%>  Specifically, OpenMPI provides the following
%>  environment variables that are defined on every MPI process:<br>
%>  <ol>
%>      <li>    ``OMPI_COMM_WORLD_SIZE``        :   The number of processes in this process ``MPI_COMM_WORLD``.<br>
%>      <li>    ``OMPI_COMM_WORLD_RANK``        :   The MPI rank of this process in ``MPI_COMM_WORLD``.<br>
%>      <li>    ``OMPI_COMM_WORLD_LOCAL_SIZE``  :   The number of ranks from this job that are running on this node.<br>
%>      <li>    ``OMPI_COMM_WORLD_LOCAL_RANK``  :   The relative rank of this process on this node within its job.<br>
%>                                                  For example, if four processes in a job share a node,
%>                                                  they will each be given a local rank ranging from ``0`` to ``3``.
%>      <li>    ``OMPI_UNIVERSE_SIZE``          :   The number of process slots allocated to this job.<br>
%>                                                  Note that this may be different than the number of processes in the job.<br>
%>      <li>    ``OMPI_COMM_WORLD_NODE_RANK``   :   The relative rank of this process on this node looking across all jobs.<br>
%>  </ol>
%>
%>  \return
%>  ``nproc``   :   The output scalar MATLAB non-negative whole number containing
%>                  the number of MPI images in the current invocation of the
%>                  OpenMPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible OpenMPI library or failed to infer the image count.<br>
%>  ``rankp1``  :   The output scalar MATLAB non-negative whole number containing
%>                  the image ID (e.g., MPI rank + 1) of the current MPI image in
%>                  the current invocation of the OpenMPI ``mpiexec`` binary MPI launcher.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible MPI library or failed to infer the image count.<br>
%>                  Note that the image ID always starts at ``1``, unlike the MPI rank.<br>
%>                  The argument ``rankp1`` stands for ``rank + 1``.<br>
%>                  An output value of ``0`` implies the algorithm failed to detect
%>                  a ParaMonte-compatible OpenMPI library or failed to infer the image count.<br>
%>
%>  \interface{ompi}
%>  \code{.m}
%>
%>      [nproc, rankp1] = pm.lib.mpi.runtime.ompi()
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
%>  \example{ompi}
%>  \include{lineno} example/lib/mpi/runtime/ompi/main.m
%>  \output{ompi}
%>  \include{lineno} example/lib/mpi/runtime/ompi/main.out.m
%>
%>  \final{ompi}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, 12:10 AM Wednesday, November 13, 2024, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function [nproc, rankp1] = ompi()

    nproc = getenv("OMPI_COMM_WORLD_SIZE");
    if ~isempty(nproc)
        nproc = str2double(nproc);
    else
        nproc = 0;
    end

    rankp1 = getenv("OMPI_COMM_WORLD_RANK");
    if ~isempty(rankp1)
        rankp1 = str2double(rankp1) + 1;
    else
        rankp1 = 0;
    end

end