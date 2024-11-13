%>  \brief
%>  Return the ID (MPI rank + 1) of the current MPI image (process),
%>  **starting from the number one**.<br>
%>
%>  \details
%>  Beware that this routine returns ``1``
%>  for the MPI rank of the first process which is ``0``.<br>
%>  This is the reason for the naming used for this function ``rankp1``
%>  standing for ``rank + 1``, consistent with the parallel ParaMonte
%>  library conventions across all supported programming languages.<br>
%>
%>  This function can be used to automate the inclusion or exclusion of
%>  single-process activities within the ParaMonte library by creating a fence
%>  around the activities (e.g., visualization, or postprocessing tasks) to
%>  be done only by a specific (normally the first) parallel process.<br>
%>  See the documentation of [pm.sampling.Paradram](@ref Paradram) for example usage.<br>
%>
%>  \param[in]  vendor  :   See the corresponding argument of [pm.lib.mpi.runtime.detect](@ref detect).<br>
%>                          (**optional**. The default is set by [pm.lib.mpi.runtime.detect](@ref detect).)
%>
%>  \return
%>  ``value``           :   The output MATLAB scalar positive whole-number containing the
%>                          image (process) ID of the MPI processes launched by the ``mpiexec``.<br>
%>                          If no MPI parallelism is involved with current application session,
%>                          of if an error occurs while inferring the MPI image ID,
%>                          the returned value is always ``1``.<br>
%>
%>  \interface{rankp1}
%>  \code{.m}
%>
%>      value = pm.lib.mpi.runtime.rankp1()
%>      value = pm.lib.mpi.runtime.rankp1([])
%>      value = pm.lib.mpi.runtime.rankp1(vendor)
%>
%>  \endcode
%>
%>  \example{rankp1}
%>  \include{lineno} example/lib/mpi/runtime/rankp1/main.m
%>  \output{rankp1}
%>  \include{lineno} example/lib/mpi/runtime/rankp1/main.out.m
%>
%>  \final{rankp1}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function value = rankp1(vendor)
    if  nargin < 1
        vendor = [];
    end
    [~, ~, value] = pm.lib.mpi.runtime.detect(vendor);
    value = max(1, value);
end