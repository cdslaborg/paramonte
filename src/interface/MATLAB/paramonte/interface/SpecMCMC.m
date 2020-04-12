classdef SpecMCMC < SpecBase

    properties
        chainSize                            = [];
        scaleFactor                          = [];
        startPointVec                        = [];
        proposalModel                        = [];
        proposalStartCovMat                  = [];
        proposalStartCorMat                  = [];
        proposalStartStdVec                  = [];
        sampleRefinementCount                = [];
        sampleRefinementMethod               = [];
        randomStartPointRequested            = [];
        randomStartPointDomainLowerLimitVec  = [];
        randomStartPointDomainUpperLimitVec  = [];
    end

end % classdef SpecMCMC
