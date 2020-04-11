classdef ParaDRAM < ParaMCMC_class

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

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    function self = ParaDRAM()

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

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    runSampler(self, ndim, getLogFunc)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    runKernel(self, getLogFunc)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    reportProgress(self)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end % methods (dynamic)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function burninLoc = getBurninLoc(lenLogFunc, refLogFunc, LogFunc)
            negLogIncidenceProb = log(lenLogFunc);
            burninLoc = 0;
            while true
                burninLoc = burninLoc + 1;
                if burninLoc < lenLogFunc && (refLogFunc - LogFunc(burninLoc)) > negLogIncidenceProb, continue; end
                break;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

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