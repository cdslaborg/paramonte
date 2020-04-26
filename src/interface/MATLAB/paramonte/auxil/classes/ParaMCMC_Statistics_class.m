classdef ParaMCMC_Statistics_class < ParaMonteStatistics_class

    properties
        Chain       = Moment_class()
        BurninLoc   = ParaMCMC_Chain_class()
        LogFuncMode = ParaMCMC_LogFuncMode_class()
    end

end