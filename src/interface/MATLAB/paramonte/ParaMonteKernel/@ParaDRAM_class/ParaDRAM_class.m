classdef ParaDRAM_class < ParaMCMC_class
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
%   Minimal serial example
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
%   Parallel simulations
%   --------------------
%
%   To run ParaDRAM sampler in parallel mode visit: cdslab.org/pm
%   You can also use the following commands on the Python command-line,
%
%       import paramonte as pm
%       pm.verify()
%
%   to obtain specific information on how to run a parallel simulation,
%   in particular, in relation to your current installation of ParaMonte.
%   In general, for parallel simulations:
%
%   0.  ensure you need and will get a speedup by running the ParaDRAM sampler in parallel.
%       Typically, if a single evaluation of the objective function takes much longer than
%       a few milliseconds, your simulation may then benefit from the parallel run.
%
%   1.  ensure you have an MPI library installed, preferably, Intel MPI.
%       An MPI library should be automatically installed on your system with ParaMonte.
%       If needed, you can download the Intel MPI library from their website and install it.
%
%   2.  set the input keyword argument mpiEnabled in runSampler() to True (the default is False),
%
%   3.  before running the parallel simulation, in particular, on Windows systems, you may
%       need to define the necessary MPI environmental variables on your system.
%       To get information on how to define the variables, use paramonte function
%       verify() as described in the above.
%
%   4.  call your main Python code from a Python-aware mpiexec-aware command-line via,
%
%           mpi_launcher -n num_process python name_of_yor_python_code.py
%
%       where,
%
%       1.  mpi_launcher is the name of the MPI launcher
%           from the MPI library that you have installed.
%           For example, the Intel MPI library launcher is named mpiexec,
%           also recognized by Microsoft and OpenMPI.
%
%       2.  num_process: replace this with the number of
%           cores on which you want to run the program.
%           Do not assign more processes than the number of
%           available physical cores on your device/cluster.
%           Assigning more cores than physically available on
%           your system will only slow down your simulation.
%
%   Minimal parallel example
%   ------------------------
%
%   Here is a Python script main.py for a parallel ParaDRAM simulation:
%
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      from paramonte import ParaDRAM
%      pmpd = ParaDRAM()
%      def getLogFunc(Point):
%          # return the natural logarithm of the ndim-dimensional
%          # non-normalized Standard Gaussian density function
%          return -sum(Point**2)
%      pmpd.runSampler ( ndim = 1
%                      , getLogFunc = getLogFunc
%                      , mpiEnabled = True
%                      )
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   where,
%
%       ndim
%           represents the number of dimensions of the domain
%           of the user's objective function getLogFunc() and,
%
%       getLogFunc()
%           represents the user's objective function that is to be sampled.
%           This function must take a single input argument of type numpy
%           float64 array of length ndim and must return the natural
%           logarithm of the objective function.
%
%
%       mpiEnabled
%           is a logical (boolean) indicator that, if True, will
%           cause the ParaDRAM simulation to run in parallel
%           on the requested number of processors.
%
%   Once the script is saved in the file main.py, open a Python-aware and
%   MPI-aware command prompt to run the simulation by the MPI launcher,
%
%       mpiexec -n 3 python main.py
%
%   This will execute the Python script main.py on three processes (images).
%   Keep in mind that on Windows systems you may need to define MPI environmental
%   variables before a parallel simulation, as descibed in the above.
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

    properties (Constant, Access = protected)
        SUB_CLASS2_NAME     = "@ParaDRAM_class"
    end

    properties (Access = public, Hidden)
        SpecDRAM            = []
        Proposal            = []
    end

    properties (Access = public)
        spec                = []
        Stats               = ParaDRAM_Statistics_class()
        Chain               = []
        RefinedChain        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = ParaDRAMKernel()

            addpath(genpath(pwd));

            self.RefinedChain   = ParaDRAMRefinedChain_class();

            self.spec           = spec();

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

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        runKernel(self, getLogFunc)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        reportProgress(self)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (dynamic)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static)

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

    methods (Hidden)    % These methods have been implemented to override the default 'handle' class methods, so that they won't pop-up after pressing 'Tab' button

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function addlistener    (self)  end
        function delete         (self)  end
        function findobj        (self)  end
        function findprop       (self)  end
        function valid          (self)  end
        function listener       (self)  end
        function notify         (self)  end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
