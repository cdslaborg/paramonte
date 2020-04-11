classdef ParaDRAM_NumFunCall_class < ParaMonteNumFunCall_class

    properties
        acceptedRejectedDelayed         = 0     % accepted + rejected function calls, including the delayed rejections
        acceptedRejectedDelayedUnused   = 0     % by all processes, used or unused
    end

end