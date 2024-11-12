%>  \brief
%>  Return the runtime number of MPI cores
%>  with which the `mpiexec` launcher may have been invoked.<br>
%>
%>  \details
%>  Otherwise, return `1` if no use of `mpiexec` launcher is detected
%>  or it is invoked with only ``1`` MPI process.<br>
%>
%>  An output value of ``1`` can be used as an indication of
%>  the ``mpiexec`` launcher in launching the ParaMonte library.<br>
%>
%>  \warning
%>  This routine can lead to a full MATLAB session crash if the
%>  required MPI library dependencies are not detected on the system.<br>
%>  This issue severely limits the utility of this routine.<br>
%>
%>  \param[in]  vendor  :   See the corresponding argument of [pm.lib.mpi.detect](@ref detect).<br>
%>                          (**optional**. The default is set by [pm.lib.mpi.detect](@ref detect).)
%>
%>  \return
%>  ``value``           :   The output MATLAB scalar positive whole number, containing the number
%>                          of MPI processes launched by the ``mpiexec`` command or ``1``
%>                          if no ``mpiexec`` invocation has occurred or ``mpiexec``
%>                          has been launched with only one process or the rank
%>                          inference fails.<br>
%>
%>  \interface{nproc}
%>  \code{.m}
%>
%>      value = pm.lib.mpi.nproc()
%>      value = pm.lib.mpi.nproc([])
%>      value = pm.lib.mpi.nproc(vendor)
%>
%>  \endcode
%>
%>  \final{nproc}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function value = nproc(vendor)
    if  nargin < 1
        vendor = [];
    end
    [~, ~, value, ~] = pm.lib.mpi.detect(vendor);
    value = max(1, value);
    %if  nargin < 1
    %    vendorList = [];
    %else
    %    vendorList = string(vendor);
    %end
    %if  pm.array.len(vendorList) == 0
    %    vendorList = pm.lib.mpi.choices();
    %end
    %value = 1;
    %mexname = "pm_parallelism";
    %for iname = 1 : length(vendorList)
    %    mexdirs = pm.lib.path.mexdir(mexname, vendorList(iname));
    %    for ipath = 1 : length(mexdirs)
    %        %mexdirs{ipath}
    %        matpath = path;
    %        pm.lib.path.clean();
    %        path(matpath, mexdirs{ipath});
    %        %persistent mh;
    %        %if ~isa(mh, 'matlab.mex.MexHost') || ~isvalid(mh)
    %        %    mh = mexhost;
    %        %end
    %        %temp = feval(mh, mexname, 'getImageCountMPI');
    %        %value = max(value, temp);
    %        mexcall = string(mexname + "('getImageCountMPI');");
    %        try
    %            temp = eval(mexcall);
    %            value = max(value, temp);
    %        catch
    %            continue;
    %        end
    %        path(matpath);
    %    end
    %end
end