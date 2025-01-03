%>  \brief
%>  Perform the basic runtime checks for the sampler and return nothing.
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
%>  \interface{run}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>      sampleList = sampler.run(getLogFunc, ndim);
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
                self.spec.parallelismNumThread = -maxNumCompThreads();
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
                    + pm.io.tab() + "self.spec.parallelismNumThread = " + string(self.parallelismNumThread) + newline ...
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
            %self.nml = [convertStringsToChars(self.getInputFile())];
            self.nml = "&" + self.method + " " + self.spec.getEntriesNML(ndim) + " " + "/" + newline;
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

    pm.lib.path.clean();
    self.matpath = path;
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
    %%%% Set up the MEX file to call.
    %%%%

    if  self.partype == "openmp"
        mexcall = string(self.mexname + "(convertStringsToChars(self.method), @getLogFuncConcurrent, ndim, convertStringsToChars(self.nml))");
        if ~self.silent
            delete(gcp("nocreate"));
            % The following works only in MATLAB 2022b and beyond.
            if  pm.matlab.release() < "2022b"
                pool = parpool("threads");
                maxNumCompThreads(abs(self.spec.parallelismNumThread));
            else
                pool = parpool("threads", abs(self.spec.parallelismNumThread));
            end
        else
            evalc('delete(gcp("nocreate"))');
            if  pm.matlab.release() < "2022b"
                evalc('pool = parpool("threads")');
                evalc('maxNumCompThreads(abs(self.spec.parallelismNumThread))');
            else
                evalc('pool = parpool("threads", abs(self.spec.parallelismNumThread))');
            end
        end
    else
        mexcall = string(self.mexname + "(convertStringsToChars(self.method), @getLogFuncWrapped, ndim, convertStringsToChars(self.nml))");
        %getLogFuncSpec = functions(getLogFunc);
    end

    %%%%
    %%%% Define the ``getLogFunc`` wrapper function for serial/MPI sampling.
    %%%%

    function logFunc = getLogFuncWrapped(state)
        logFunc = getLogFunc(state);
    end

    %%%%
    %%%% Define the ``getLogFunc`` wrapper function for OpenMP sampling.
    %%%%

    %getLogFuncConst = parallel.pool.Constant(@(state) getLogFunc(state));
    function [logFunc, avgTimePerFunCall, avgCommPerFunCall] = getLogFuncConcurrent(state)
        %state = parallel.pool.Constant(state(2 : end, :));
        avgCommPerFunCall = tic();
        avgTimePerFunCall = 0;
        njob = size(state, 2);

        % %getState = arrayfun(@(ijob) state(:, ijob), 1 : njob, 'UniformOutput', 0);
        % fout(1 : njob) = parallel.FevalFuture;
        % start = zeros(njob, 'uint64');
        % delta = zeros(njob);
        % for ijob = 1 : njob
        %     %slice = getState(ijob);
        %     start(ijob) = tic();
        %     %fout(ijob) = parfeval(pool, getLogFunc, 1, slice{1});
        %     %fout(ijob) = parfeval(getLogFunc, 1, state(:, ijob));
        %     fout(ijob) = parfeval(pool, getLogFunc, 1, state(:, ijob));
        %     delta(ijob) = toc(start(ijob));
        % end
        % avgTimePerFunCall = sum(delta);
        % logFunc = fetchOutputs(fout);

        %spmd
        %    indx = spmdIndex;
        %    avgTimePerFunCall = tic();
        %    logFunc = getLogFunc(state(:, indx));
        %    avgTimePerFunCall = toc(avgTimePerFunCall);
        %end
        %logFunc = [logFunc{:}];
        %avgTimePerFunCall = mean([avgTimePerFunCall{:}]);

        %%pause

        % \devnote
        % `parfor` is slightly slower than `parfeval` (abour 20 ms for density function costing .1 ms).
        % However, the slight penalty diminishes for costly problems and furthermore, `parfor` allows
        % a reliable method of timing the density function, unlike `parfeval`.

        logFunc = zeros(njob, 1);
        %getSlice = arrayfun(@(ijob) state(:, ijob), 1 : njob, 'UniformOutput', 0);
        parfor ijob = 1 : njob
            start = tic();
            %logFunc(ijob) = getLogFunc(getSlice{ijob});
            %logFunc(ijob) = getLogFunc(state(:, ijob));
            logFunc(ijob) = feval(getLogFunc, state(:, ijob));
            avgTimePerFunCall = avgTimePerFunCall + toc(start);
        end

        avgTimePerFunCall = avgTimePerFunCall / njob;
        avgCommPerFunCall = toc(avgCommPerFunCall) - avgTimePerFunCall;
    end

    %%%%
    %%%% Call the MEX sampler.
    %%%%

    try
        eval(mexcall);
        self.finalize();
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
        self.finalize();
        msg = string(me.identifier) + " : " + string(me.message) + newline;
        if ismac && strcmpi(me.identifier, 'MATLAB:mex:ErrInvalidMEXFile')
            msg = msg ...
                + "This error is most likely due to the ""System Integrity Protection""" + newline ...
                + "(SIP) of your macOS interfering with the ParaMonte MATLAB MEX files." + newline ...
                + "You can follow the guidelines in the documentation to resolve this error:" + newline ...
                + newline ...
                + pm.io.tab() + pm.web.href(self.weblinks.generic.docs.url + "/troubleshooting/macos-developer-cannot-be-verified/") + newline ...
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
                + pm.io.tab() + self.spec.outputFileName + "*_report.txt" + newline ...
                + newline ...
                ;
        end
        error(msg);
    end

end