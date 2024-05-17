function count = nproc(config)
    %
    %   Return the runtime number of MPI cores
    %   with which the `mpiexec` launcher may have been invoked.
    %   Otherwise, return `1` if no use of `mpiexec` launcher is detected
    %   or it is invoked with only ``1`` MPI process.
    %
    %   An output value of ``1`` can be used as an indication of
    %   the ``mpiexec`` launcher in launching the ParaMonte library.
    %
    %   \warning
    %
    %       This routine can lead to a full MATLAB session crash
    %       if the required MPI library dependencies are not detected
    %       on the system. This issue severely limits the utility of this routine.
    %
    %   Parameters
    %   ----------
    %
    %       config
    %
    %           The input scalar (or array of) MATLAB string(s)
    %           containing the ParaMonte-preferred MPI library vendor/name
    %           or other configurations as used in naming the ParaMonte MATLAB shared libraries.
    %           This argument is passed directly to the corresponding argument of ``pm.lib.path.mexdir``.
    %           Possible values of outmost interest to MPI applications are:
    %
    %               -   ``Intel`` or ``impi``, representing the Intel MPI library.
    %               -   ``MPICH`` or `mmpi``, representing the MPICH MPI library.
    %               -   ``OpenMPI`` or ``ompi``, representing the OpenMPI library.
    %
    %           All values are case-insensitive.
    %           (**optional**. If missing, all possibilities are considered and the largest inferred ``count`` is returned.)
    %
    %   Returns
    %   -------
    %
    %       count
    %
    %           The output MATLAB scalar integer containing the number of
    %           MPI processes launched by the ``mpiexec`` command or ``1``
    %           if no ``mpiexec`` invocation has occurred or the routine
    %           fails to load any MPI-enabled ParaMonte library.
    %
    %   Interface
    %   ---------
    %
    %       count = pm.lib.mpi.nproc()
    %       count = pm.lib.mpi.nproc([])
    %       count = pm.lib.mpi.nproc(config)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if nargin < 1
        configList = [];
    else
        configList = string(config);
    end
    if  pm.array.len(configList) == 0
        configList = pm.lib.mpiChoices();
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