%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonteSampler - This is base class for the ParaMonte sampler routines.
%   
%   Parameters
%   ----------
%
%       platform
%
%           A MATLAB struct containing the platform logical values (isMacOS, isWin32, isLinux).
%           Note that ll class attributes can be set after an instance
%           is returned by the constructor.
%
%   Attributes
%   ----------
%
%       See below for information on the attributes (properties).  
%
%   Methods
%   -------
%
%       See below for information on the methods.  
%   
%   Returns
%   -------
%
%       Object of class ParaMonteSampler
%
%   Naming conventions
%   ------------------
%
%       camelCase naming style is used throughout the entire ParaMonte library, across
%       all programming languages: C/Fortran/Julia/MATLAB/Python
%
%       All simulation specifications start with a lowercase letter, including
%       scalar/vector/matrix int, float, string, or boolean variables.
%
%       The name of any variable that represents a vector of values is suffixed with "Vec",
%       for example: startPointVec, domainLowerLimitVec, ...
%
%       The name of any variable that represents a matrix of values is suffixed with "Mat",
%       for example: proposalStartCorMat, ...
%
%       The name of any variable that represents a list of varying-size values is suffixed
%       with "List", for example: variableNameList, ...
%
%       All functions or class methods begin with a lowercase verb.
%
%       Significant attempt has been made to end all boolean variables with a passive verb,
%       such that the full variable name virtually forms an English-language statement
%       that should be either true or false, set by the user.
%
%   Tips
%   ----
%
%       When running ParaMonte samplers, in particular on multiple cores in parallel,
%       it would be best to close any such aggressive software/applications as
%       Dropbox, ZoneAlarm, ... that can interfere with your ParaMonte
%       simulation output files, potentially causing the sampler to
%       crash before successful completion of the simulation.
%       These situations should however happen only scarcely.
%
%       On Windows systems, when restarting an old interrupted ParaMonteSampler simulation,
%       ensure your MATLAB session is also restarted before the simulation
%       restart. This is needed as Windows sometimes locks the access to
%       all or some of the simulation output files.
%
%       To unset an already-set input simulation specification, simply set the
%       simulation attribute to None or re-instantiate the sampler object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef ParaMonteSampler < dynamicprops

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        %
        %       buildMode
        %           optional string argument with default value "release".
        %           possible choices are:
        %               "debug"
        %                   to be used for identifying sources of bug
        %                   and causes of code crash.
        %               "release"
        %                   to be used in all other normal scenarios
        %                   for maximum runtime efficiency.
        buildMode = "release";
        %
        %       mpiEnabled
        %           optional logical (boolean) indicator which is ``false`` by default.
        %           If it is set to ``true``, it will cause the ParaMonte simulation
        %           to run in parallel on the requested number of processors.
        %           See the class documentation guidelines in the above for 
        %           information on how to run a simulation in parallel.
        mpiEnabled = false;
        %
        %       reportEnabled
        %           optional logical (boolean) indicator which is ``true`` by default.
        %           If it is set to ``false``, it will silence all output postprocessing
        %           messages. If ``mpiEnabled=true``, ``reportEnabled=false`` automatically.
        reportEnabled = true;
        %
        %       inputFile
        %           optional string input representing the path to
        %           an external input namelist of simulation specifications.
        %           USE THIS OPTIONAL ARGUMENT WITH CAUTION AND
        %           ONLY IF YOU KNOW WHAT YOU ARE DOING.
        %           ==================================================================
        %           WARNING: Specifying an input file will cause the ParaMonte sampler 
        %           to ignore all other simulation specifications set by the user via 
        %           sampler instance's `spec`-component attributes.
        %           ==================================================================
        inputFile = [];
        %
        %       spec
        %           A MATLAB structure containing all simulation specifications.
        %           All simulation attributes are by default set to appropriate 
        %           values at runtime. To override the default simulation 
        %           specifications, set the `spec` attributes to some 
        %           desired values of your choice.
        spec = [];
    end

    properties (Access = public, Hidden)
        Err = Err_class();
        methodName = "";
        objectName = [];
        platform = [];
        libName = [];
        website = [];
        method = [];
        ndim = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = ParaMonteSampler(platform,website)
            %filePath = mfilename("fullpath"); addpath(genpath(filePath),"-begin");
            self.method = struct();
            self.method.isParaDRAM = false;
            self.method.isParaNest = false;
            self.method.isParaTemp = false;
            self.website = website;
            self.platform = platform;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        chainList = readChain(self,varargin)
        sampleList = readSample(self,varargin)
        reportList = readReport(self,varargin)
        progressList = readProgress(self,varargin)
        restartList = readRestart(self,varargin)
        runSampler(self,ndim,getLogFunc,varargin)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)
        % These methods have been implemented to override the default 'handle' class methods, 
        % so that they won't pop-up after pressing 'Tab' button.
        function addlistener    (); end
        function delete         (); end
        function findobj        (); end
        function findprop       (); end
        function valid          (); end
        function listener       (); end
        function notify         (); end
        function name = getMyName(self); name = inputname(1); end
        delFile(self, file, desc)
        [filePathList, iswebfile] = getFilePathList(self,file,fileType)
        result = genOutputFileName(self)
        namelist = getInputFile(self)
        outputList = readOutput(self,file,delimiter,fileType)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef ParaMonteSampler
