classdef SpecMCMC < SpecBase

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

end % classdef SpecMCMC
