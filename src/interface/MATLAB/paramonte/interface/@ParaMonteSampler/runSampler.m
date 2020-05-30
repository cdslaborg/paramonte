%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see,
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms,
%   if you use any parts of this library for any purposes,
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   runSampler(ndim,getLogFunc)
%
%   Run the kernel sampler and return nothing.
%
%   Parameters
%   ----------
%
%       ndim
%           integer representing the number of dimensions of the
%           domain of the user's objective function getLogFunc().
%           It must be a positive integer.
%
%       getLogFunc()
%           represents the user's objective function to be sampled,
%           which must take a single input argument of type numpy
%           float64 array of length ndim and must return the
%           natural logarithm of the objective function.
%
%   Returns
%   -------
%
%       None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function runSampler(self,ndim,getLogFunc,varargin)
    if nargin~=3
        self.Err.msg    = "The method " + self.objectName + ".runSampler(ndim,getLogFunc) takes only two input arguments:" + newline + newline ...
                        + "          ndim:  the number of dimensions of the domain of the " + newline ...
                        + "                 mathematical objective function to be sampled," + newline ...
                        + "    getLogFunc:  a MATLAB handle to the MATLAB function implementing " + newline ...
                        + "                 the mathematical objective function to be sampled,";
        self.Err.abort();
    end

    self.ndim = ndim;

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

    if ~isa(getLogFunc,"function_handle")
        %if exist("getLogFunc") && ( isa(YourVar,'function_handle')
        %    if isequal(getLogFunc,
        %    end
        %end
    %else
        self.Err.msg    = "The input argument getLogFunc must be a callable function. " + newline ...
                        + "It represents the user's objective function to be sampled, " + newline ...
                        + "which must take a single input argument of type numpy " + newline ...
                        + "float64 array of length ndim and must return the " + newline ...
                        + "natural logarithm of the objective function.";
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errorOccurred = ~isa(self.buildMode,"string");
    cstype = [];
    if ~errorOccurred
        buildModeSplitList = strsplit(lower(self.buildMode),"-");
        if strcmp(buildModeSplitList(1),"release") || strcmp(buildModeSplitList(1),"testing") || strcmp(buildModeSplitList(1),"debug")
            if length(buildModeSplitList) > 1
                cstype = buildModeSplitList(2);
                errorOccurred = ~( strcmp(cstype,"gnu") || strcmp(cstype,"intel") );
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

    if self.mpiEnabled
        parallelism = "_mpi";
    else
        parallelism = "";
        self.Err.msg    = "Running the ParaDRAM sampler in serial mode..." + newline ...
                        + "To run the ParaDRAM sampler in parallel mode visit: cdslab.org/pm";
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
    if ~isempty(cstype)
        pmcsList = [cstype];
        for pmcs = pmcsListRef
            if ~strcmp(pmcs,cstype)
                pmcsList = [ pmcsList , [pmcs] ];
            end
        end
    end

    libFound = false;
    if self.platform.isWin32; libNameSuffix = "_windows_" + getArch() + "_mt"; end
    if self.platform.isMacOS; libNameSuffix = "_darwin_" + getArch() + "_mt"; end
    if self.platform.isLinux; libNameSuffix = "_linux_" + getArch() + "_mt"; end
    for buildMode = buildModeList
        for pmcs = pmcsList
            libName = "libparamonte_dynamic_heap_" + buildMode + "_" + pmcs + parallelism + libNameSuffix;
            if exist(libName,'file')==3; libFound = true; break; end
        end
        if libFound; break; end
    end
    self.libName = libName;

    if ~libFound
        self.Err.msg    = "Exhausted all possible ParaMonte dynamic library search" + newline ...
                        + "names but could not find any compatible library." + newline ...
                        + "It appears your ParaMonte Python interface is missing" + newline ...
                        + "the dynamic libraries. Please report this issue at:" + newline + newline ...
                        + "    <a href=""https://github.com/cdslaborg/paramonte/issues"">https://github.com/cdslaborg/paramonte/issues</a>" + newline + newline ...
                        + "Visit <a href=""https://www.cdslab.org/paramonte/"">https://www.cdslab.org/paramonte/</a> for instructions" + newline ...
                        + "to build ParaMonte object files on your system."; %...+ newline ...
                        %...+ buildInstructionNote % xxxxxxxxxxxx...
        self.Err.abort();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function logFunc = getLogFuncNested(point)
        logFunc = getLogFunc(point);
    end
    getLogFuncSpec = functions(getLogFunc);

    try
        iscmd = isdeployed() || batchStartupOptionUsed; % batchStartupOptionUsed is introduced in R2019a and not supported in older versions of MATLAB
    catch
        iscmd = isdeployed();
    end
    %if ~(self.mpiEnabled || iscmd)
    %    self.Err.msg = "check the opened terminal for simulation progress and report.";
    %    self.Err.note();
    %end

    if strcmp(getLogFuncSpec.type,"simple") && strcmp(getLogFuncSpec.function,"getLogFunc")
        expression = string(self.libName + "(iscmd,ndim,inputFile)");
    else
        expression = string(self.libName + "(iscmd,ndim,inputFile,@getLogFuncNested)");
    end

    isGNU = contains(self.libName,"gnu");
    if isGNU
        setenv('GFORTRAN_STDIN_UNIT' , '5') 
        setenv('GFORTRAN_STDOUT_UNIT', '6') 
        setenv('GFORTRAN_STDERR_UNIT', '0')
    end
    munlock(self.libName)
    eval("clear "+self.libName);
    eval(expression);
    munlock(self.libName)
    eval("clear "+self.libName);
    if isGNU
        setenv('GFORTRAN_STDIN_UNIT' , '-1') 
        setenv('GFORTRAN_STDOUT_UNIT', '-1') 
        setenv('GFORTRAN_STDERR_UNIT', '-1')
    end

    %try
    %    eval(expression);
    %catch
    %    if self.mpiEnabled
    %        reportFileSuffix = "_process_*_report.txt";
    %    else
    %        reportFileSuffix = "_process_1_report.txt";
    %    end
    %    self.Err.msg    = "The " + self.methodName + " simulation failed. Please see the contents of the simulation's output report file(s) for potential diagnostic messages:" + newline + newline ...
    %                    + "    " + strrep(self.spec.outputFileName+reportFileSuffix,'\','\\') ...
    %                    ;
    %    self.Err.abort();
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~self.mpiEnabled
        self.Err.msg    = "To read the generated output files sample or chain files, try the following:" + newline + newline ...
                        + "    " + self.objectName + ".readSample()      %% to read the final i.i.d. sample from the output sample file." + newline + newline ...
                        + "    " + self.objectName + ".readChain()       %% to read the uniquely-accepted points from the output chain file." + newline + newline ...
                        + "    " + self.objectName + ".readMarkovChain() %% to read the Markov Chain. NOT recommended for extremely-large chains." + newline + newline + newline ...
                        + "For more information and examples on the usage, visit:" + newline + newline ...
                        + "    <a href=""https://www.cdslab.org/paramonte/"">https://www.cdslab.org/paramonte/</a>";
        self.Err.note();
    end

end