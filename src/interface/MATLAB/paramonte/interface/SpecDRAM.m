classdef SpecDRAM < SpecMCMC

    properties
        adaptiveUpdateCount                  = [];
        adaptiveUpdatePeriod                 = [];
        greedyAdaptationCount                = [];
        delayedRejectionCount                = [];
        burninAdaptationMeasure              = [];
        delayedRejectionScaleFactorVec       = [];
    end

end % classdef SpecDRAM
