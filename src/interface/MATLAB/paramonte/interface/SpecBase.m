classdef SpecBase < handle

    properties
        sampleSize                          = []
        randomSeed                          = []
        description                         = []
        outputFileName                      = []
        outputDelimiter                     = []
        chainFileFormat                     = []
        variableNameList                    = []
        restartFileFormat                   = []
        outputColumnWidth                   = []
        outputRealPrecision                 = []
        silentModeRequested                 = []
        domainLowerLimitVec                 = []
        domainUpperLimitVec                 = []
        parallelizationModel                = []
        progressReportPeriod                = []
        targetAcceptanceRate                = []
        maxNumDomainCheckToWarn             = []
        maxNumDomainCheckToStop             = []
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        % These methods have been implemented to override the default 'handle' class methods, 
        % so that they won't pop-up after pressing 'Tab' button.

        function addlistener    (self)  end
        function delete         (self)  end
        function findobj        (self)  end
        function findprop       (self)  end
        function valid          (self)  end
        function listener       (self)  end
        function notify         (self)  end

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecBase
