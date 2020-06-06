%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it 
%   under the terms of the GNU Lesser General Public License as published 
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see, 
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms, 
%   if you use any parts of this library for any purposes, 
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef Global_class < handle

    properties
        timeElapsedUntilLastReportInSeconds     = []
        SumAccRateSinceStart                    = []    
        sumAccRateLastReport                    = []
        inverseProgressReportPeriod             = []
        numFunCallAcceptedRejectedLastReport    = []
        numFunCallAccepted_dummy                = []
        co_proposalFound_samplerUpdateOccurred  = []
        comv_chol                               = []
        comv_covMat                             = []
        mc_restartFileUnit                      = []
        mv_sampleSizeOld_save                   = []
        mv_logSqrtDetOld_save                   = []
        mv_adaptiveScaleFactorSq_save           = []
        mv_MeanOld_save                         = []
        mc_Image                                = []     
        mc_ndim                                 = []     
        mc_logFileUnit                          = []     
        mc_scalingRequested                     = []     
        mc_defaultScaleFactorSq                 = []     
        mc_DelayedRejectionCount                = []     
        mc_MaxNumDomainCheckToWarn              = []     
        mc_MaxNumDomainCheckToStop              = []     
        mc_delayedRejectionRequested            = []     
        mc_ndimInverse                          = []     
        mc_targetAcceptanceRate                 = []     
        mc_DelayedRejectionScaleFactorVec       = []     
        mc_DomainLowerLimitVec                  = []     
        mc_DomainUpperLimitVec                  = []     
        mc_MaxNumDomainCheckToWarnMsg           = []     
        mc_MaxNumDomainCheckToStopMsg           = []     
        mc_negativeHellingerDistSqMsg           = []     
        mc_restartFileFormat                    = []     
        mc_methodBrand                          = []     
        mc_methodName                           = []     
        mc_isNormal                             = []
        mv_Err                                  = []
    end

end