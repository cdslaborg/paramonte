%>  \brief
%>  Validate the simulation settings, execute the selected native ParaMonte sampler,
%>  and finalize its run-scoped MATLAB resources. This method returns nothing.
%>
%>  \param[inout]   self        :   The input/output parent object of class [pm.sampling.Sampler](@ref Sampler)
%>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
%>  \param[in]      getLogFunc  :   The input MATLAB function handle or anonymous (lambda) function
%>                                  containing the implementation of the objective function to be sampled.<br>
%>                                  This user-specified function must have the following interface,<br>
%>                                  \code{.m}
%>                                      function logFunc = getLogFunc(state)
%>                                          ...
%>                                      end<br>
%>                                  \endcode
%>                                  where,<br>
%>                                  <ol>
%>                                      <li>    the input argument ``state`` is a vector of type MATLAB ``double``
%>                                              of size ``ndim`` representing a single point from within the ``ndim``
%>                                              dimensional domain of the mathematical object function to be explored.<br>
%>                                      <li>    the output argument `logFunc` is a scalar of the same type as the
%>                                              input ``state`` containing the natural logarithm of the objective
%>                                              function at the specified input ``state`` within its domain.<br>
%>                                  </ol>
%>  \param[in]      ndim        :   The input scalar positive-valued whole-number representing the number of dimensions
%>                                  of the domain of the user-specified objective function in the input ``getLogFunc()``.<br>
%>
%>  \details
%>  MPI takes precedence over thread parallelism. If MPI is not selected,
%>  ``spec.parallelismNumThread = []`` selects serial execution; a positive
%>  integer selects thread parallelism, and zero uses ``maxNumCompThreads()``.
%>  Thread-parallel objective dispatch is selected automatically: MATLAB R2025a
%>  and later use persistent queue workers; earlier releases with thread-pool
%>  support use timed ``parfor`` on a reused thread pool. No additional sampler
%>  property or simulation specification is required.<br>
%>
%>  The native OpenMP callback passes an ``ndim``-by-``njob`` matrix of
%>  candidate states to the dispatcher. Both backends return an ``njob``-by-1
%>  vector in original column order, the arithmetic mean of the individual
%>  evaluation times (seconds), and the nonnegative difference between batch
%>  wall time and that mean (seconds). The latter is not pure communication time.<br>
%>
%>  In persistent mode, one resident ``parfeval`` task per worker remains active
%>  until the simulation ends. A private closure owned by this run manages
%>  worker state; neither the closure nor its handles are exposed as a public
%>  sampler property. The objective handle and its captured data are
%>  passed at startup, not on each callback. The callback is synchronous: the
%>  native sampler resumes only after all proposals in the batch have completed.<br>
%>
%>  \note
%>  Both paths require Parallel Computing Toolbox and a thread-capable objective.
%>  A compatible thread pool is reused and retained after the run, even if created
%>  by the sampler. R2022b and later require the requested worker count; older
%>  releases use the existing/default count with a warning if it differs.
%>  Persistent dispatch requires an idle pool. An existing incompatible pool
%>  is never deleted or resized automatically. Only the resident futures
%>  created for this run are stopped. The cleanup guard removes
%>  the temporary MEX path and restores saved GNU Fortran environment settings.
%>  No output argument is returned; use the sampler's output-reading methods
%>  to access the generated files.<br>
%>
%>  \warning
%>  The objective must be safe for concurrent execution on MATLAB thread workers.
%>  Do not rely on execution order, worker-local mutable state, or random-stream
%>  placement; do not submit nested work to the occupied pool. Worker errors abort
%>  the batch and are propagated to the caller. Cleanup on interruption is best
%>  effort for uninterruptible native or objective code.<br>
%>
%>  \note
%>  If ``input`` supplies an external namelist, its native thread configuration
%>  must agree with ``spec.parallelismNumThread``. The internal dispatcher
%>  selection does not change the native namelist or callback interface.<br>
%>  \interface{run}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      sampler.spec.parallelismNumThread = 4;
%>      sampler.run(getLogFunc, ndim);
%>
%>  \endcode
%>
%>  \final{run}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:38 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function run(self, getLogFunc, ndim)

    %%%%
    %%%% Sanitize ``sampler.silent``.
    %%%%

    if ~pm.introspection.verified(self.silent, "logical", 1)
        help("pm.sampling.Sampler.silent");
        disp("self.silent =");
        disp(self.silent);
        error   ( newline ...
                + "The sampler attribute ``silent`` must be a scalar of type ``logical``." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end

    %%%%
    %%%% Sanitize parallelism method to set reporting permission.
    %%%%

    % global mpiname;
    % ismember('mpiname', who('global'));
    % MPI enabled by a global definition of ``mpiname``.
    % if  pm.array.len(self.mpiname) == 0 && ~isempty(mpiname) && pm.introspection.verified(mpiname, "string", 1)
    %     self.mpiname = mpiname;
    % end

    if ~pm.introspection.verified(self.mpiname, "string", 1)

        %%%%
        %%%% Sanitize ``mpiname``.
        %%%%

        help("pm.sampling.Sampler.mpiname");
        disp("mpiname =");
        disp(self.mpiname);
        error   ( newline ...
                + "The sampler attribute ``mpiname`` must be a scalar of type ``char`` or ``string``." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );

    elseif pm.array.len(self.mpiname) == 0

        %%%%
        %%%% Detect potential MPI launcher.
        %%%%

        [mpiname, nproc, ~] = pm.lib.mpi.runtime.detect();
        if  pm.array.len(mpiname) > 0 && nproc > 1
            self.mpiname = mpiname;
        end

    end

    %%%%
    %%%% First detect potential use of MPI, then check for thread parallelism.
    %%%%

    self.partype = "serial";
    numThread = [];
    if  0 < pm.array.len(self.mpiname) % MPI enabled.

        self.silent = true; % Otherwise, we keep the default value of self.silent.
        self.partype = string(pm.lib.mpi.name(self.mpiname));

    elseif ~isempty(self.spec.parallelismNumThread)

        %%%%
        %%%% Sanitize ``self.spec.parallelismNumThread``.
        %%%%

        % The following separate conditions are crucial to remain separate.
        failed = ~pm.introspection.verified(self.spec.parallelismNumThread, "integer", 1);
        if ~failed
            failed = self.spec.parallelismNumThread < 0;
        end
        if ~failed
            failed = ~pm.matlab.has.parallel;
        end
        if ~failed
            self.partype = "openmp";
            if  self.spec.parallelismNumThread == 0
                %%%% Negative threads count indicates to the ParaMonte sampler that the user set `parallelismNumThread` to zero.
                numThread = maxNumCompThreads();
            else
                numThread = double(self.spec.parallelismNumThread);
            end
        else
            help("pm.sampling.Sampler.run");
            error   ( newline ...
                    + "The simulation specification ``parallelismNumThread`` must be" + newline ...
                    + "a scalar non-negative integer or whole-number representing the number" + newline ...
                    + "of processor threads to use for a thread-parallel simulation." + newline ...
                    + "Specifying ``0`` as the value of ``parallelismNumThread`` will" + newline ...
                    + "lead to using all available CPU processes for thread-parallelism." + newline ...
                    + "Ideally, a number less than the maximum available number" + newline ...
                    + "of threads should be specified for ``parallelismNumThread`` as using" + newline ...
                    + "all available threads exclusively for a simulation will" + newline ...
                    + "slow down all open system applications including MATLAB." + newline ...
                    + "Beware that specifying this option requires the MATLAB parallel toolbox." + newline ...
                    + "If missing or specified as empty ``[]``, the simulation will run in serial mode." + newline ...
                    + "You have specified:" + newline ...
                    + newline ...
                    + pm.io.tab() + "self.spec.parallelismNumThread = " + string(self.spec.parallelismNumThread) + newline ...
                    + newline ...
                    + "Does your MATLAB have Parallel Computing Toolbox?" + newline ...
                    + newline ...
                    + pm.io.tab() + "pm.matlab.has.parallel() = " + string(pm.matlab.has.parallel()) + newline ...
                    + newline ...
                    + "For more information, see the documentation displayed above." + newline ...
                    + newline ...
                    );
        end

    end

    %%%%
    %%%% Sanitize ``getLogFunc``.
    %%%%

    if ~isa(getLogFunc, "function_handle")
        help("pm.sampling.Sampler.run");
        error   ( newline ...
                + "The input argument getLogFunc must be a callable function." + newline ...
                + "It represents the user's objective function to be sampled," + newline ...
                + "which must take a single input argument of type numpy" + newline ...
                + "float64 array of length ndim and must return the" + newline ...
                + "natural logarithm of the objective function." + newline ...
                + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end

    %%%%
    %%%% Sanitize ``ndim``.
    %%%%

    failed = ~pm.introspection.verified(ndim, "integer", 1);
    if ~failed
        failed = ndim < 1;
    end
    if ~failed
        self.ndim = ndim;
    else
        help("pm.sampling.Sampler.run");
        disp("ndim = ");
        disp(ndim);
        error   ( newline ...
                + "The input argument ``ndim`` must be a scalar positive integer (or whole-number)." + newline ...
                + "The input argument ``ndim`` represents the number of dimensions of the domain" + newline ...
                + "of the user-specified objective function ``getLogFunc()``." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    end

    %%%%
    %%%% Sanitize ``input`` specifications/file string.
    %%%%

    if ~pm.introspection.verified(self.input, "string", 1)
        help("pm.sampling.Sampler.input");
        disp("self.input = ");
        disp(self.input);
        error   ( newline ...
                + "The specified sampler attribute ``input`` must be scalar MATLAB string." + newline ...
                + "For more information, see the documentation displayed above." + newline ...
                + newline ...
                );
    else
        self.input = string(self.input);
        if  self.input ~= ""
            if isfile(self.input)
                self.nml = pm.sys.path.abs(self.input, "lean");
            else
                self.nml = self.input; % Could still be a valid namelist.
            end
            if ~self.silent
                warning ( newline ...
                        + "User-specified input namelist file detected: " + newline ...
                        + newline ...
                        + pm.io.tab() + """" + self.input + """" + newline ...
                        + newline ...
                        + "All simulation specifications will be read from the input file." + newline ...
                        + "All simulation specifications in the ``spec`` component of the sampler object will be ignored." + newline ...
                        + newline ...
                        );
            end
        else
            % Preserve the native negative-count sentinel without permanently
            % changing the user-facing specification from zero to a negative value.
            if self.partype == "openmp" && self.spec.parallelismNumThread == 0
                entries = getEntriesWithThreadSentinel(self.spec, ndim, numThread);
            else
                entries = self.spec.getEntriesNML(ndim);
            end
            self.nml = "&" + self.method + " " + entries + " " + "/" + newline;
        end
    end

    %%%%
    %%%% Sanitize ``checked``.
    %%%%

    if ~isempty(self.checked)
        if ~pm.introspection.verified(self.checked, "logical", 1)
            help("pm.sampling.Sampler.checked");
            disp("self.checked =");
            disp(self.checked);
            error   ( newline ...
                    + "The specified sampler attribute ``checked`` must be a scalar MATLAB logical." + newline ...
                    + "For more information, see the documentation displayed above." + newline ...
                    + newline ...
                    );
        elseif ~self.checked
            chktypes = "nocheck";
        else
            chktypes = "checked";
        end
    else
        chktypes = ["nocheck", "checked"];
    end

    %%%%
    %%%% Setup the ParaMonter sampler library name.
    %%%%

    libspecs = [pm.os.namel(), pm.sys.arch(), self.libtype, self.memtype, self.partype];
    mexdirs = pm.lib.path.mexdir(self.mexname, libspecs);

    % We will choose the checking mode, compiler suite,
    % and build mode based on the availability below.

    if  self.bldtype == ""
        bldtypes = pm.lib.bldtypes();
    else
        bldtypes = lower(self.bldtype);
    end
    if  self.clstype == ""
        clstypes = pm.lib.clstypes();
    else
        clstypes = lower(self.clstype);
    end
    failed = isempty(mexdirs);
    if ~failed
        failed = true;
        for ichk = 1 : length(chktypes)
            chktype = filesep + chktypes(ichk);
            for icls = 1 : length(clstypes)
                clstype = filesep + clstypes(icls);
                for ibld = 1 : length(bldtypes)
                    bldtype = filesep + bldtypes(ibld) + filesep;
                    for imex = 1 : length(mexdirs)
                        if  contains(mexdirs(imex), clstype) ...
                        &&  contains(mexdirs(imex), bldtype) ...
                        &&  contains(mexdirs(imex), chktype)
                            mexdir = mexdirs(imex);
                            if  self.bldtype == ""
                                self.bldtype = bldtypes(ibld);
                            end
                            if  self.clstype == ""
                                self.clstype = clstypes(icls);
                            end
                            failed = false;
                            break;
                        end
                    end
                    if ~failed
                        break;
                    end
                end
                if ~failed
                    break;
                end
            end
            if ~failed
                break;
            end
        end
    end

    if  failed
        help("pm.sampling.Sampler");
        disp("libspecs");
        disp( libspecs );
        disp("bldtypes");
        disp( bldtypes );
        disp("clstypes");
        disp( clstypes );
        error   ( newline ...
                + "There are no MEX libraries associated with the configurations displayed above:" + newline ...
                + "Either the user has compromised internal structure of the ParaMonte library" + newline ...
                + "or the user has tempered with hidden attributes of the ParaMonte sampler." + newline ...
                + "If you believe neither is the case, please report this error at:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(self.weblinks.github.issues.url) + newline ...
                + newline ...
                + "for a quick resolution." + newline ...
                + newline ...
                );
    end

    %%%%
    %%%% Add the identified MEX path to the MATLAB path list, only temporarily.
    %%%%

    dispatcher = [];
    if self.partype == "openmp"
        dispatcher = getObjectiveDispatcher(@serveObjectiveJobs, @evaluateObjectiveParfor);
    end

    pm.lib.path.clean();
    self.matpath = path;
    originalEnv = {getenv('GFORTRAN_STDIN_UNIT'), ...
                   getenv('GFORTRAN_STDOUT_UNIT'), ...
                   getenv('GFORTRAN_STDERR_UNIT')};
    runGuard = onCleanup(@() self.finalize(dispatcher, originalEnv));
    path(self.matpath, mexdir);
    munlock(self.mexname);
    clear(self.mexname);

    if ~self.silent && ~(pm.matlab.iscmd || pm.os.is.win)
        disp( newline ...
            + "Also check the shell terminal (from which you opened MATLAB)" + newline ...
            + "for potential realtime simulation progress and report." + newline ...
            + newline ...
            );
    end

    if  self.clstype == "gnu"
        setenv('GFORTRAN_STDIN_UNIT' , '5');
        setenv('GFORTRAN_STDOUT_UNIT', '6');
        setenv('GFORTRAN_STDERR_UNIT', '0');
    end

    %%%%
    %%%% Initialize the objective dispatcher and call the MEX sampler.
    %%%%

    try
        if self.partype == "openmp"
            dispatcher.start(numThread, getLogFunc, self.silent);
        end
        invokeMex(self, getLogFunc, ndim, dispatcher);
        % Stop the resident workers before restoring paths or printing the report.
        clear runGuard;
        if ~self.silent
            disp( newline ...
                + self.getppm() + newline ...
                + "For more information and examples on the usage, visit:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(self.weblinks.docs.generic.url) + newline ...
                + newline ...
                );
        end
    catch me
        clear runGuard;
        msg = string(me.identifier) + " : " + string(me.message) + newline;
        if ismac && strcmpi(me.identifier, 'MATLAB:mex:ErrInvalidMEXFile')
            msg = msg ...
                + "This error is most likely due to the ""System Integrity Protection""" + newline ...
                + "(SIP) of your macOS interfering with the ParaMonte MATLAB MEX files." + newline ...
                + "You can follow the guidelines in the documentation to resolve this error:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(self.weblinks.docs.generic.url + "/troubleshooting/macos-developer-cannot-be-verified/") + newline ...
                + newline ...
                + "If the problem persists even after following the guidelines" + newline ...
                + "in the above page, please report this issue to the developers at:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(self.weblinks.github.issues.url) ...
                + newline ...
                ;
        else
            msg = msg ...
                + "The " + self.method + " simulation failed. See the error message above." + newline ...
                + "Also check the contents of the generated output '*_report.txt' files" + newline ...
                + "if any such files were generated before the simulation crash:" + newline ...
                + newline ...
                + pm.io.tab() + reportPrefix(self.spec.outputFileName) + "*_report.txt" + newline ...
                + newline ...
                ;
        end
        failure = MException('pm:sampling:Sampler:RunFailed', '%s', char(msg));
        failure = addCause(failure, me);
        throw(failure);
    end

end

function entries = getEntriesWithThreadSentinel(spec, ndim, numThread)
    originalSpec = spec.parallelismNumThread;
    specGuard = onCleanup(@() restoreThreadSpec(spec, originalSpec));
    spec.parallelismNumThread = -numThread;
    entries = spec.getEntriesNML(ndim);
end

function restoreThreadSpec(spec, originalSpec)
    spec.parallelismNumThread = originalSpec;
end

function invokeMex(self, getLogFunc, ndim, dispatcher)
    % Keep the original callback signatures, including all three OpenMP
    % outputs. This nested workspace owns no run-level cleanup objects.
    if self.partype == "openmp"
        callback = @getLogFuncConcurrent;
    else
        callback = @getLogFuncWrapped;
    end
    feval(char(self.mexname), convertStringsToChars(self.method), ...
          callback, ndim, convertStringsToChars(self.nml));

    function logFunc = getLogFuncWrapped(state)
        logFunc = getLogFunc(state);
    end

    function [logFunc, avgTimePerFunCall, avgCommPerFunCall] = getLogFuncConcurrent(state)
        [logFunc, avgTimePerFunCall, avgCommPerFunCall] = dispatcher.evaluate(state);
    end
end

function prefix = reportPrefix(outputFileName)
    if isempty(outputFileName)
        prefix = "";
    else
        prefix = string(outputFileName);
    end
end
