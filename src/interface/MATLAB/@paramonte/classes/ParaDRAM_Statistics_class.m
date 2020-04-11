classdef ParaDRAM_Statistics_class < ParaMCMC_Statistics_class

    properties
        AdaptationBurninLoc = ParaMCMC_Chain_class()
        NumFunCall          = ParaDRAM_NumFunCall_class()
    end

end