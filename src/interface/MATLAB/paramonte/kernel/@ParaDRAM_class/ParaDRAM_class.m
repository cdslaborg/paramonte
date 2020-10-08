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
%   This is the ParaDRAMKernel class for generating instances of serial and parallel
%   Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo
%   sampler of the ParaMonte library.
%
%   All ParaDRAM class attributes (input arguments to the ParaDRAM constructor)
%   are optional and all attributes can be also set after a ParaDRAM instance
%   is returned by the constructor.
%
%   Once you set the desired attributes to the desired values,
%   call the ParaDRAM sampler via the object's method runSampler().
%
%   Minimal example
%   ----------------------
%
%   Here is a MATLAB script main.m for a serial ParaDRAM simulation.
%   Copy and paste this code into your MATLAB command line to run it:
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       pmlibRootDir = '../'; % set this path to the ParaMonte library root directory
%       addpath( genpath(pmlibRootDir) );
%       %The following lambda function returns the natural logarithm of the
%       %ndim-dimensional non-normalized Standard Gaussian density function
%       getLogFunc = @(point) -0.5 * sum(point.^2);
%       pmpd = ParaDRAM();
%       pmpd.runSampler ( 4 ... number of dimensions of the objective function
%                       , getLogFunc ... the mathematical objective function
%                       );
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   where,
%
%       ndim
%           represents the number of dimensions of the domain
%           of the user's objective function getLogFunc() and,
%
%       getLogFunc()
%           represents the user's objective function to be sampled,
%           which must take a single input argument of type numpy
%           float64 array of length ndim and must return the
%           natural logarithm of the objective function.
%
%   ParaDRAM Class Attributes
%   -------------------------
%
%   (see also: cdslab.org/pm)
%
%   All input specifications (attributes) of a ParaDRAM simulation are optional.
%   However, it is recomended that you provide as much information as possible
%   about the specific ParaDRAM simulation and the objective function to be sampled
%   via ParaDRAM simulation specifications.
%
%   The ParaDRAM simulation specifications have lengthy comprehensive descriptions
%   that appear in full in the output report file of every ParaDRAM simulation.
%
%   The best way to learn about individual ParaDRAM simulation attributes
%   is to a run a minimal serial simulation with the following Python script,
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       from paramonte import ParaDRAM
%       pmpd = ParaDRAM()
%       pmpd.outputFileName = "./test"
%       def getLogFunc(Point): return -sum(Point**2)
%       pmpd.runSampler( ndim = 1, getLogFunc = getLogFunc )
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Running this code will generate a set of simulation output files (in the current
%   working directory of Python) that begin with the prefix "test_process_1". Among
%   these, the file "test_process_1_report.txt" contains the full description of all
%   input specifications of the ParaDRAM simulation as well as other information
%   about the simulation results and statistics.
%
%   Naming conventions
%   ------------------
%
%   camelCase naming style is used throughout the entire ParaMonte library, across
%   all programming languages: C/Fortran/Julia/MATLAB/Python
%
%   all simulation specifications start with a lowercase letter, including
%   scalar/vector/matrix int, float, string, or boolean variables.
%
%   The name of any variable that represents a vector of values is suffixed with "Vec",
%   for example: startPointVec, domainLowerLimitVec, ...
%
%   The name of any variable that represents a matrix of values is suffixed with "Mat",
%   for example: proposalStartCorMat, ...
%
%   The name of any variable that represents a list of varying-size values is suffixed
%   with "List", for example: variableNameList, ...
%
%   all functions or class methods begin with a lowercase verb.
%
%   significant attempt has been made to end all boolean variables with a passive verb,
%   such that the full variable name virtually forms an English-language statement
%   that should be either True or False, set by the user.
%
%   Tips
%   ----
%
%   when running ParaMonte samplers, in particular on multiple cores in parallel,
%   it would be best to close any such aggressive software/applications as
%   Dropbox, ZoneAlarm, ... that can interfere with your ParaMonte
%   simulation output files, potentially causing the sampler to
%   crash before successful completion of the simulation.
%   These situations should however happen only scarcely.
%
%   on Windows systems, when restarting an old interrupted ParaDRAM simulation,
%   ensure your Python session is also restarted before the simulation
%   restart. This is needed as Windows sometimes locks the access to
%   all or some of the simulation output files.
%
%   To unset an already-set input simulation specification, simply set the
%   simulation attribute to None or re-instantiate the object.
classdef ParaDRAM_class < ParaMCMC_class

    properties (Constant, Access = protected, Hidden)
        SUB_CLASS2_NAME     = "@ParaDRAM_class"
    end

    properties (Access = public, Hidden)
        methodName          = "MatDRAM";    % the algorithm's name
        mpiEnabled          = false;        % no parallelization implemented in MatDRAM yet
        objectName          = [];           % dynamic name of the user-defined objects
        inputFile           = [];           % no input file in MatDRAM
        SpecDRAM            = [];
        Proposal            = [];
        Stats               = ParaDRAM_Statistics_class()
        Chain               = [];
        RefinedChain        = [];
    end

    properties (Access = public)
        reportEnabled       = true;         % no verbose messaging, if false;
        spec                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = ParaDRAM_class(platform,website,version)

            addpath(genpath(pwd));

            self = self@ParaMCMC_class(platform,website);

            self.RefinedChain   = ParaDRAMRefinedChain_class();

            self.version = version;
            self.spec = spec();

            ...ParaMonte variables
            self.spec.sampleSize                            = [];
            self.spec.randomSeed                            = [];
            self.spec.description                           = [];
            self.spec.outputFileName                        = [];
            self.spec.outputDelimiter                       = [];
            self.spec.chainFileFormat                       = [];
            self.spec.variableNameList                      = [];
            self.spec.restartFileFormat                     = [];
            self.spec.outputColumnWidth                     = [];
            self.spec.overwriteRequested                    = [];
            self.spec.outputRealPrecision                   = [];
            self.spec.silentModeRequested                   = [];
            self.spec.domainLowerLimitVec                   = [];
            self.spec.domainUpperLimitVec                   = [];
            self.spec.parallelizationModel                  = [];
            self.spec.progressReportPeriod                  = [];
            self.spec.targetAcceptanceRate                  = [];
            self.spec.maxNumDomainCheckToWarn               = [];
            self.spec.maxNumDomainCheckToStop               = [];
            ...ParaMCMC variables
            self.spec.chainSize                             = [];
            self.spec.startPointVec                         = [];
            self.spec.sampleRefinementCount                 = [];
            self.spec.sampleRefinementMethod                = [];
            self.spec.randomStartPointRequested             = [];
            self.spec.randomStartPointDomainLowerLimitVec   = [];
            self.spec.randomStartPointDomainUpperLimitVec   = [];
            ...ParaDRAM variables
            self.spec.scaleFactor                           = [];
            self.spec.proposalModel                         = [];
            self.spec.proposalStartCovMat                   = [];
            self.spec.proposalStartCorMat                   = [];
            self.spec.proposalStartStdVec                   = [];
            self.spec.adaptiveUpdateCount                   = [];
            self.spec.adaptiveUpdatePeriod                  = [];
            self.spec.greedyAdaptationCount                 = [];
            self.spec.delayedRejectionCount                 = [];
            self.spec.burninAdaptationMeasure               = [];
            self.spec.delayedRejectionScaleFactorVec        = [];

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        runSampler(self, ndim, getLogFunc)
        [markovChainList] = readMarkovChain(self,varargin)
        [restartList] = readRestart(self,varargin)
        [reportList] = readReport(self,varargin)
        [sampleList] = readSample(self,varargin)
        [chainList] = readChain(self,varargin)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (dynamic)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static, Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function burninLoc = getBurninLoc(lenLogFunc, refLogFunc, LogFunc)
            negLogIncidenceProb = log(lenLogFunc);
            burninLoc = 0;
            while true
                burninLoc = burninLoc + 1;
                if burninLoc < lenLogFunc && (refLogFunc - LogFunc(burninLoc)) > negLogIncidenceProb, continue; end
                break;
            end
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Hidden)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        runKernel(self, getLogFunc)
        readTabular(self, callerName, varargin)
        reportProgress(self)
        writeOutput(self)
        writeOutput2(self)
        writeRestartFile(self)
        [filePathList, iswebfile] = getFilePathList(self,file,fileType)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        % These methods have been implemented to override the default 'handle' class methods, so that they won't pop-up after pressing 'Tab' button

        function addlistener(self)  end
        function delete     (self)  end
        function findobj    (self)  end
        function findprop   (self)  end
        function valid      (self)  end
        function listener   (self)  end
        function notify     (self)  end
        function eq         (self)  end
        function ge         (self)  end
        function gt         (self)  end
        function le         (self)  end
        function lt         (self)  end
        function ne         (self)  end

        % cannot overrise the folowing method
        % function isvalid    (self)  end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
