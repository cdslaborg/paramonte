%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        lower_comv_covMat                       = []
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