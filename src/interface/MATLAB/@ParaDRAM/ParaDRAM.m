classdef ParaDRAM < ParaMCMC_class

    properties (Constant, Access = protected)
        SUB_CLASS2_NAME     = "@ParaDRAM_class"
    end

    properties (Access = public, Hidden)
        SpecDRAM            = []
        Proposal            = []
    end

    properties (Access = public)
        specs               = []
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

        self.specs          = specs();

        ...ParaMonte variables
        self.specs.sampleSize                           = [];
        self.specs.randomSeed                           = [];
        self.specs.description                          = [];
        self.specs.outputFileName                       = [];
        self.specs.outputDelimiter                      = [];
        self.specs.chainFileFormat                      = [];
        self.specs.variableNameList                     = [];
        self.specs.restartFileFormat                    = [];
        self.specs.outputColumnWidth                    = [];
        self.specs.outputRealPrecision                  = [];
        self.specs.silentModeRequested                  = [];
        self.specs.domainLowerLimitVec                  = [];
        self.specs.domainUpperLimitVec                  = [];
        self.specs.parallelizationModel                 = [];
        self.specs.progressReportPeriod                 = [];
        self.specs.targetAcceptanceRate                 = [];
        self.specs.maxNumDomainCheckToWarn              = [];
        self.specs.maxNumDomainCheckToStop              = [];
        ...ParaMCMC variables
        self.specs.chainSize                            = [];
        self.specs.startPointVec                        = [];
        self.specs.sampleRefinementCount                = [];
        self.specs.sampleRefinementMethod               = [];
        self.specs.randomStartPointRequested            = [];
        self.specs.randomStartPointDomainLowerLimitVec  = [];
        self.specs.randomStartPointDomainUpperLimitVec  = [];
        ...ParaDRAM variables
        self.specs.scaleFactor                          = [];
        self.specs.proposalModel                        = [];
        self.specs.proposalStartCovMat                  = [];
        self.specs.proposalStartCorMat                  = [];
        self.specs.proposalStartStdVec                  = [];
        self.specs.adaptiveUpdateCount                  = [];
        self.specs.adaptiveUpdatePeriod                 = [];
        self.specs.greedyAdaptationCount                = [];
        self.specs.delayedRejectionCount                = [];
        self.specs.burninAdaptationMeasure              = [];
        self.specs.delayedRejectionScaleFactorVec       = [];

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