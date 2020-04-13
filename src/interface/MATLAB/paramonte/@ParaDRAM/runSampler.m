function runSampler(self,ndim,getLogFunc)
%   Run ParaDRAM sampler and return nothing.
%   
%   Parameters
%   ----------
%       ndim
%           integer representing the number of dimensions of the
%           domain of the user's objective function getLogFunc().
%           It must be a positive integer.
%       getLogFunc()
%           represents the user's objective function to be sampled,
%           which must take a single input argument of type numpy
%           float64 array of length ndim and must return the
%           natural logarithm of the objective function.
%   
%   Returns
%   -------
%       None

    self.objectName = inputname(1);

    self.Err.marginTop = 0;
    self.Err.marginBot = 1;
    self.Err.resetEnabled = false;

    if ~isa(ndim,"numeric") || ndim<1
        self.Err.msg    = "The input argument ndim must be a positive integer, " + newline ...
                        + "representing the number of dimensions of the domain of " + newline ...
                        + "the user's objective function getLogFunc(). " + newline ...
                        + "You have entered ndim = " + string(ndim);
        self.Err.abort();
    end

    if ~isa(getLogFunc,'function_handle')
        self.Err.msg    = "The input argument getLogFunc must be a callable function. " + newline ...
                        + "It represents the user's objective function to be sampled, " + newline ...
                        + "which must take a single input argument of type numpy " + newline ...
                        + "float64 array of length ndim and must return the " + newline ...
                        + "natural logarithm of the objective function.";
        self.Err.abort();
    end

    errorOccurred = ~isa(self.buildMode,"string");
    stype = [];
    if ~errorOccurred
        buildModeSplitList = strsplit(lower(self.buildMode),"-");
        if strcmp(buildModeSplitList(1),"release") || strcmp(buildModeSplitList(1),"testing") || strcmp(buildModeSplitList(1),"debug")
            if length(buildModeSplitList) > 1
                stype = buildModeSplitList(2);
                errorOccurred = ~( strcmp(stype,"gnu") || strcmp(stype,"intel") );
            end
        else
            errorOccurred = true;
        end
    end
    if errorOccurred
        self.Err.msg    = "The attribute " + self.objectName + ".buildMode must be of type string. " + newline ...
                        + "It is an optional string argument with default value ""release"" ." + newline ...
                        + "possible choices are: " + newline ...
                        + "    ""debug"": " + newline ...
                        + "        to be used for identifying sources of bug " + newline ...
                        + "        and causes of code crash. " + newline ...
                        + "    ""release"": " + newline ...
                        + "        to be used in all other normal scenarios " + newline ...
                        + "        for maximum runtime efficiency. " + newline ...
                        + "You have entered " + self.objectName + ".buildMode = " + string(self.buildMode);
        self.Err.abort();
    end

    if ~isa(self.mpiEnabled,"logical")
        self.Err.msg    = "The input argument " + self.objectName + ".mpiEnabled must be of type bool. " + newline ...
                        + "It is an optional logical (boolean) indicator which is False by default. " + newline ...
                        + "If it is set to True, it will cause the ParaDRAM simulation " + newline ...
                        + "to run in parallel on the requested number of processors. " + newline ...
                        + "See ParaDRAM class information on how to run a simulation in parallel. " + newline ...
                        + "You have entered " + self.objectName + ".mpiEnabled = " + string(self.mpiEnabled);
        self.Err.abort();
    end

    if ~isempty(self.inputFile) && ~isa(self.inputFile,"string")
        self.Err.msg    = "The input argument " + self.objectName + ".inputFile must be of type string. " + newline ...
                        + "It is an optional string input representing the path to " + newline ...
                        + "an external input namelist of simulation specifications. " + newline ...
                        + "USE THIS OPTIONAL ARGUMENT WITH CAUTION AND " + newline ...
                        + "ONLY IF YOU KNOW WHAT YOU ARE DOING. " + newline ...
                        + "Specifying this option will cause ParaDRAM to ignore " + newline ...
                        + "all other paraDRAM simulation specifications set by " + newline ...
                        + "the user via ParaDRAM instance attributes. " + newline ...
                        + "You have entered " + self.objectName + ".inputFile = """ + string(self.inputFile) + """.";
        self.Err.warn();
    end

    inputFile = [convertStringsToChars(self.getInputFile())];
    %inputFile = self.getInputFile()

    isWin32 = ispc;
    isMacOS = ismac;
    isLinux = isunix;

    if self.mpiEnabled
        parallelism = "_mpi";
    else
        parallelism = "";
        self.Err.msg    = "Running the ParaDRAM sampler in serial mode..." + newline ...
                        + "To run the ParaDRAM sampler in parallel mode visit: cdslab.org/pm" + newline ...
                        + "check the opened terminal for simulation progress and report.";
        self.Err.note();
    end
    buildModeListRef = ["release","testing","debug"];
    buildModeList = buildModeListRef;
    if ~strcmp(self.buildMode,"release")
        for buildMode = buildModeListRef
            if ~strcmp(buildMode,self.buildMode)
                buildModeList = [ buildModeList , [buildMode] ];
            end
        end
    end

    pmcsListRef = ["intel","gnu"];
    pmcsList = pmcsListRef;
    if ~isempty(stype)
        pmcsList = [stype];
        for pmcs = pmcsListRef
            if ~strcmp(pmcs,stype)
                pmcsList = [ pmcsList , [pmcs] ];
            end
        end
    end

    libFound = false;
    if isWin32; libnameSuffix = "_windows_" + getArch() + "_mt"; end
    if isMacOS; libnameSuffix = "_darwin_" + getArch() + "_mt"; end
    if isLinux; libnameSuffix = "_linux_" + getArch() + "_mt"; end
    for buildMode = buildModeList
        for pmcs = pmcsList
            libname = "libparamonte_dynamic_heap_" + buildMode + "_" + pmcs + "_m" + parallelism + libnameSuffix;
            if exist(libname,'file')==3; libFound = true; break; end
        end
        if libFound; break; end
    end

    if ~libFound
        self.Err.msg    = "Exhausted all possible ParaMonte dynamic library search" + newline ...
                        + "names but could not find any compatible library." + newline ...
                        + "It appears your ParaMonte Python interface is missing" + newline ...
                        + "the dynamic libraries. Please report this issue at:" + newline + newline ...
                        + "    https://github.com/cdslaborg/paramonte/issues" + newline + newline ...
                        + "Visit https://www.cdslab.org/pm for instructions" + newline ...
                        + "to build ParaMonte object files on your system."; ...+ newline ...
                        ...+ buildInstructionNote % xxxxxxxxxxxx...
        self.Err.abort();
    end

    expression = string(libname + "(ndim,inputFile,isdeployed())");
    eval(expression);

    if ~self.mpiEnabled
        self.Err.msg    = "To read the generated output files sample or chain files, try the following:" + newline + newline ...
                        + "    " + self.objectName + ".readSample()      # to read the final i.i.d. sample from the output sample file." + newline ...
                        + "    " + self.objectName + ".readChain()       # to read the uniquely-accepted points from the output chain file." + newline ...
                        + "    " + self.objectName + ".readMarkovChain() # to read the Markov Chain. NOT recommended for extremely-large chains." + newline + newline ...
                        + "For more information and examples on the usage, visit:" + newline + newline ...
                        + "    https://www.cdslab.org/paramonte/";
        self.Err.note();
    end

end