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
%>  \param[in]  config  :   The input scalar (or array of) MATLAB string(s)
%>                          containing the ParaMonte-preferred MPI library vendor/name
%>                          or other configurations as used in naming the ParaMonte MATLAB shared libraries.<br>
%>                          This argument is passed directly to the corresponding argument of ``pm.lib.path.mexdir``.<br>
%>                          Possible values of outmost interest to MPI applications are:<br>
%>                          <ol>
%>                              <li>    ``Intel`` or ``impi``, representing the Intel MPI library.<br>
%>                              <li>    ``MPICH`` or `mmpi``, representing the MPICH MPI library.<br>
%>                              <li>    ``OpenMPI`` or ``ompi``, representing the OpenMPI library.<br>
%>                          </ol>
%>                          Note that **all values are case-insensitive**.<br>
%>                          (**optional**. If missing, all possibilities are considered and the largest inferred ``count`` is returned.)
%>
%>  \return
%>  `count`             :   The output MATLAB scalar integer containing the number of
%>                          MPI processes launched by the ``mpiexec`` command or ``1``
%>                          if no ``mpiexec`` invocation has occurred or the routine
%>                          fails to load any MPI-enabled ParaMonte library.<br>
%>
%>  \interface{nproc}
%>  \code{.m}
%>
%>      count = pm.lib.mpi.nproc()
%>      count = pm.lib.mpi.nproc([])
%>      count = pm.lib.mpi.nproc(config)
%>
%>  \endcode
%>
%>  \final{nproc}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:14 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function count = nproc(config)
    if  nargin < 1
        configList = [];
    else
        configList = string(config);
    end
    if  pm.array.len(configList) == 0
        configList = pm.lib.mpi.choices();
    end
    count = 1;
    mexname = "pm_parallelism";
    for iname = 1 : length(configList)
        mexdirs = pm.lib.path.mexdir(mexname, configList(iname));
        for ipath = 1 : length(mexdirs)
            %mexdirs{ipath}
            matpath = path;
            pm.lib.path.clean();
            path(matpath, mexdirs{ipath});
            %persistent mh;
            %if ~isa(mh, 'matlab.mex.MexHost') || ~isvalid(mh)
            %    mh = mexhost;
            %end
            %temp = feval(mh, mexname, 'getImageCountMPI');
            %count = max(count, temp);
            mexcall = string(mexname + "('getImageCountMPI');");
            try
                temp = eval(mexcall);
                count = max(count, temp);
            catch
                continue;
            end
            path(matpath);
        end
    end
end