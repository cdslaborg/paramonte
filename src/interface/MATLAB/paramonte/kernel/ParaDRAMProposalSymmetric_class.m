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

classdef ParaDRAMProposalSymmetric_class < handle

    properties (Constant)
        SUB_CLASS_NAME      = "@ParaDRAMProposalSymmetric_class"
    end

    properties
       %AccRate             = struct('sumUpToLastUpdate', [], 'target', [])
        Err                 = Err_class()
        Global              = Global_class()
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = ParaDRAMProposalSymmetric_class(ndim, PD)

            FUNCTION_NAME = self.SUB_CLASS_NAME + "@ParaDRAMProposalSymmetric_class()";

            %***********************************************************************************************************************
            % Global variables
            %***********************************************************************************************************************

%            global  comv_chol                           ...
%                    comv_covMat                         ...
%                    mc_Image                            ...
%                    mc_ndim                             ...
%                    mc_logFileUnit                      ...
%                    mc_restartFileUnit                  ...
%                    mc_scalingRequested                 ...
%                    mc_defaultScaleFactorSq             ...
%                    mc_DelayedRejectionCount            ...
%                    mc_MaxNumDomainCheckToWarn          ...
%                    mc_MaxNumDomainCheckToStop          ...
%                    mc_delayedRejectionRequested        ...
%                    mc_ndimInverse                      ...
%                    mc_targetAcceptanceRate             ...
%                    mc_DelayedRejectionScaleFactorVec   ...
%                    mc_DomainLowerLimitVec              ...
%                    mc_DomainUpperLimitVec              ...
%                    mc_MaxNumDomainCheckToWarnMsg       ...
%                    mc_MaxNumDomainCheckToStopMsg       ...
%                    mc_negativeAdaptationMeasureMsg       ...
%                    mc_restartFileFormat                ...
%                    mc_methodBrand                      ...
%                    mc_methodName                       ...
%                    mc_isNormal                         ...
%                    ...the following had to be defined globally for the sake of restart file generation
%                    mv_MeanOld_save                     ...
%                    mv_logSqrtDetOld_save               ...
%                    mv_adaptiveScaleFactorSq_save       ...
%                    mv_sampleSizeOld_save               ...
%                    mv_Err

            %***********************************************************************************************************************
            % setup sampler update global save variables
            %***********************************************************************************************************************

            self.Global.mv_MeanOld_save(1:ndim)         = 0; % Constants.NULL_RK;
            self.Global.mv_logSqrtDetOld_save           = 0; % Constants.NULL_RK;
            self.Global.mv_sampleSizeOld_save           = 0;
            self.Global.mv_adaptiveScaleFactorSq_save   = 1;
            self.Global.mv_Err                          = Err_class();

            %***********************************************************************************************************************
            % setup general ParaDRAMProposalSymmetric specifications
            %***********************************************************************************************************************

            self.Global.mc_ndim                             = ndim;
            self.Global.mc_DomainLowerLimitVec              = PD.SpecBase.domainLowerLimitVec.Val;
            self.Global.mc_DomainUpperLimitVec              = PD.SpecBase.domainUpperLimitVec.Val;
            self.Global.mc_DelayedRejectionScaleFactorVec   = PD.SpecDRAM.delayedRejectionScaleFactorVec.Val;
            self.Global.mc_isNormal                         = PD.SpecDRAM.proposalModel.isNormal;
            self.Global.mc_Image                            = PD.Image;
            self.Global.mc_methodName                       = PD.name;
            self.Global.mc_methodBrand                      = PD.brand;
            self.Global.mc_logFileUnit                      = PD.LogFile.unit;
            self.Global.mc_restartFileUnit                  = PD.RestartFile.unit;
            self.Global.mc_restartFileFormat                = PD.RestartFile.format;
            % self.Global.mc_restartFileFormat                = PD.SpecBase.restartFileFormat.val;
            self.Global.mc_defaultScaleFactorSq             = PD.SpecDRAM.scaleFactor.val^2;
           %self.AccRate.sumUpToLastUpdate      = 0;
            self.Global.mc_MaxNumDomainCheckToWarn          = PD.SpecBase.maxNumDomainCheckToWarn.val;
            self.Global.mc_MaxNumDomainCheckToStop          = PD.SpecBase.maxNumDomainCheckToStop.val;
            self.Global.mc_DelayedRejectionCount            = PD.SpecDRAM.delayedRejectionCount.val;
            self.Global.mc_delayedRejectionRequested        = self.Global.mc_DelayedRejectionCount > 0;
            self.Global.mc_scalingRequested                 = PD.SpecBase.targetAcceptanceRate.scalingRequested;

            if self.Global.mc_scalingRequested
                self.Global.mc_targetAcceptanceRate         = PD.SpecBase.targetAcceptanceRate.val;
            else
                self.Global.mc_targetAcceptanceRate         = 1.0;    % 0.234
            end

            self.Global.mc_ndimInverse                      = 1.0 / ndim;

            %***********************************************************************************************************************
            % setup ProposalSymmetric specifications
            %***********************************************************************************************************************

            % setup covariance matrix

            self.Global.comv_covMat = zeros(ndim, ndim, self.Global.mc_DelayedRejectionCount + 1); % xxx: define new variable 
            self.Global.lower_comv_covMat   = zeros(1, ndim*(ndim-1)/2 + ndim);
            self.Global.comv_chol   = zeros(ndim, ndim, self.Global.mc_DelayedRejectionCount + 1);

            if PD.SpecDRAM.proposalStartCovMat.isPresent
                self.Global.comv_covMat(:, :, 1)  = PD.SpecDRAM.proposalStartCovMat.Val;
            else
                self.Global.comv_covMat(:, :, 1)  = getCovMatFromCorMat ( ndim                                  ...
                                                            , PD.SpecDRAM.proposalStartStdVec.Val   ...
                                                            , PD.SpecDRAM.proposalStartCorMat.Val   ...
                                                            ) ;
            end

            % Now scale the covariance matrix

            self.Global.comv_covMat(:,:,1) = self.Global.comv_covMat(:,:,1) * self.Global.mc_defaultScaleFactorSq;

            % Now get the Cholesky Factor of the Covariance Matrix. Lower comv_CholDiagLower will be the CholFac

            self.Global.comv_chol(:,:,1)    = chol(self.Global.comv_covMat(:, :, 1), 'lower'); % xxx check what happens if chol fails

            % if comv_chol(1, 1, 1) < 0
                % mv_Err.msg    = mc_Image.name + FUNCTION_NAME + ": Singular input covariance matrix by user was detected. " ...
                                % + "This is strange." + newline + "Covariance matrix lower triangle:";
                % for j = 1 : ndim
                    % for i = 1 : j
                        % mv_Err.msg = mv_Err.msg + newline + mat2str(comv_covMat(1:i, j, 1));
                    % end
                % end
                % mv_Err.abort(mc_methodBrand, newline, mc_logFileUnit);
            % end

            % Scale the higher-stage delayed-rejection Cholesky Lower matrices

            if self.Global.mc_delayedRejectionRequested, self.updateDelRejCholDiagLower(); end

            % This will be used for Domain boundary checking during the simulation

            self.Global.mc_negativeAdaptationMeasureMsg   = self.SUB_CLASS_NAME + "@doAdaptation(): Non-positive adaptation measure detected, "           ...
                                            + "possibly due to round-off error: ";

            self.Global.mc_MaxNumDomainCheckToWarnMsg   = self.SUB_CLASS_NAME + "@getNew(): " + num2str(self.Global.mc_MaxNumDomainCheckToWarn) + " proposals were "...
                                            + "drawn out of the objective function's Domain without any acceptance.";

            self.Global.mc_MaxNumDomainCheckToStopMsg   = self.SUB_CLASS_NAME + "@getNew(): " + num2str(self.Global.mc_MaxNumDomainCheckToStop)                     ...
                                            + " proposals were drawn out of the objective function's Domain. As per the value set "         ...
                                            + "for the simulation specification variable 'maxNumDomainCheckToStop', " + self.Global.mc_methodName       ...
                                            + " will abort now.";

        end % ParaDRAMProposalSymmetric_class

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function StateNew = getNew  ( self          ...
                                    , nd            ...
                                    , counterDRS    ...
                                    , StateOld      ...
                                    )
            StateOld = transpose(StateOld);
    
            % FUNCTION_NAME = ParaDRAMProposalSymmetric_class.SUB_CLASS_NAME + "@getNew()";

%            global  comv_chol                       ...
%                    ... comv_covMat                     ...
%                    mc_logFileUnit                  ...
%                    mc_MaxNumDomainCheckToWarn      ...
%                    mc_MaxNumDomainCheckToStop      ...
%                    mc_DomainLowerLimitVec          ...
%                    mc_DomainUpperLimitVec          ...
%                    mc_MaxNumDomainCheckToWarnMsg   ...
%                    mc_MaxNumDomainCheckToStopMsg   ...
%                    mc_methodBrand                  ...
%                    mc_isNormal                     ...
%                    mv_Err
                    
            self.Global.mv_Err.prefix       = self.Global.mc_methodBrand;
            self.Global.mv_Err.outputUnit   = self.Global.mc_logFileUnit;

            domainCheckCounter = 0;

            if self.Global.mc_isNormal

                while true  % loopBoundaryCheckNormal: Check for the support Region consistency:

                    % StateNew = mvnrnd(StateOld, comv_chol(:,:,counterDRS+1));

                    StateNew = getRandMVN   ( nd                            ...
                                            , StateOld                      ...
                                            , self.Global.comv_chol (:,:,counterDRS+1)  ...
                                            ) ;

                    %StateNew = StateOld' + comv_chol (:,:,counterDRS+1)*randn(nd,1);
                    %StateNew = StateNew';

                    if any( StateNew(:) < self.Global.mc_DomainLowerLimitVec(:) ) || any( StateNew(:) > self.Global.mc_DomainUpperLimitVec(:) )
                        domainCheckCounter = domainCheckCounter + 1;
                        if domainCheckCounter == self.Global.mc_MaxNumDomainCheckToWarn
                            self.Global.mv_Err.msg      = self.Global.mc_MaxNumDomainCheckToWarnMsg;
                            self.Global.mv_Err.warn();
                        end
                        if domainCheckCounter == self.Global.mc_MaxNumDomainCheckToStop
                            self.Global.mv_Err.msg      = self.Global.mc_MaxNumDomainCheckToStopMsg;
                            self.Global.mv_Err.abort();
                        end
                        continue
                    end
                    break
                end % loopBoundaryCheckNormal

            else % if (mc_isUniform)

                while true  % loopBoundaryCheckUniform

                    StateNew = getRandMVU   ( nd                            ...
                                            , StateOld                      ...
                                            , self.Global.comv_chol (:,:,counterDRS+1)  ...
                                            ) ;

                    if any(StateNew <= self.Global.mc_DomainLowerLimitVec) || any(StateNew >= self.Global.mc_DomainUpperLimitVec)
                        domainCheckCounter  = domainCheckCounter + 1;
                        if domainCheckCounter == self.Global.mc_MaxNumDomainCheckToWarn
                            self.Global.mv_Err.msg      = self.Global.mc_MaxNumDomainCheckToWarnMsg;
                            self.Global.mv_Err.warn();
                        end
                        if domainCheckCounter == self.Global.mc_MaxNumDomainCheckToStop
                            self.Global.mv_Err.msg      = self.Global.mc_MaxNumDomainCheckToStopMsg;
                            self.Global.mv_Err.abort();
                        end
                        continue
                    end
                    break

                end % loopBoundaryCheckUniform

            end % if

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************
    
        % ATTN: This routine may need further correction for the delayed rejection method
        function adaptationMeasure =  doAutoTune( self              ...
                                                , adaptationMeasure   ...
                                                , AutoTuneScaleSq   ...
                                                )
%            global  comv_chol               ...
%                    comv_covMat             ...
%                    mv_logSqrtDetOld_save   ...
%                    mv_Err
                                            
            FUNCTION_NAME           = SUB_CLASS_NAME + "@doAutoTune()";

            self.Global.mv_Err.prefix           = self.Global.mc_methodBrand;
            self.Global.mv_Err.outputUnit       = self.Global.mc_logFileUnit;

            CovMatUpperOld          = self.Global.comv_covMat(:,:,1);
            self.Global.mv_logSqrtDetOld_save   = sum(log(diag(self.Global.comv_chol(:,:,1))));

            if AutoTuneScaleSq(1) == 0
                self.Global.comv_covMat(1,1,1)  = 0.25 * self.Global.comv_covMat(1,1,1);
                self.Global.comv_chol(1,1,1)    = sqrt(self.Global.comv_covMat(1,1,1));
            else
                self.Global.comv_covMat(1,1,1)  = AutoTuneScaleSq(1);
                self.Global.comv_chol(1,1,1)    = sqrt(AutoTuneScaleSq(1));
            end

            % compute the adaptationMeasure (adaptivity)
            logSqrtDetNew       = sum(log(diag(self.Global.comv_chol(:,:,1))));
            CovMatUpperCurrent  = 0.5 * (self.Global.comv_covMat(:,:,1) + CovMatUpperOld);
            [logSqrtDetSum, singularityOccurred] = getLogSqrtDetPosDefMat(1, CovMatUpperCurrent);
            if singularityOccurred
                self.Global.mv_Err.msg  = FUNCTION_NAME                                                                             ...
                            + ": Error occurred while computing the Cholesky factorization of "                         ...
                            + "a matrix needed for the computation of the proposal distribution's adaptation measure. " ...
                            + "Such error is highly unusual, and requires an in depth investigation of the case. "      ...
                            + "Restarting the simulation might resolve the error."                                      ...
                            ;
                self.Global.mv_Err.abort();
                return
            end
            adaptationMeasure = 1 - exp(0.5 * (self.Global.mv_logSqrtDetOld_save + logSqrtDetNew) - logSqrtDetSum);

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [samplerUpdateSucceeded, adaptationMeasure]  = doAdaptation( self                      ...
                                                                            , nd                        ...
                                                                            , chainSize                 ...
                                                                            , Chain                     ...
                                                                            , ChainWeight               ...
                                                                            , samplerUpdateIsGreedy     ...
                                                                            , meanAccRateSinceStart     ...
                                                                            , adaptationMeasure         ...
                                                                            )
            FUNCTION_NAME = ParaDRAMProposalSymmetric_class.SUB_CLASS_NAME + "@doAdaptation()";

%            global  comv_chol                           ...
%                    comv_covMat                         ...
%                    mc_logFileUnit                      ...
%                    mc_scalingRequested                 ...
%                    mc_defaultScaleFactorSq             ...
%                    mc_delayedRejectionRequested        ...
%                    mc_ndimInverse                      ...
%                    mc_targetAcceptanceRate             ...
%                    mc_negativeadaptationMeasureMsg       ...
%                    mc_methodBrand                      ...
%                    ...the following had to be defined globally for the sake of restart file generation
%                    mv_MeanOld_save                     ...
%                    mv_logSqrtDetOld_save               ...
%                    mv_adaptiveScaleFactorSq_save       ...
%                    mv_sampleSizeOld_save               ...
%                    mv_Err

            scalingNeeded           = false;
            ... sampleSizeOld           = mv_sampleSizeOld_save;  % this is kept only for restoration of mv_sampleSizeOld_save, if needed.
            self.Global.mv_logSqrtDetOld_save   = sum(log(diag(self.Global.comv_chol(:,:,1))));

            % First if there are less than nd+1 points for new covariance computation, then just scale the covariance and return

            if chainSize > nd   % blockSufficientSampleSizeCheck

                % get the new sample's upper covariance matrix and Mean

                if samplerUpdateIsGreedy
                    sampleSizeCurrent = chainSize;
                    [CovMatUpperCurrent, MeanCurrent] = getSamCovUpperMeanTrans (sampleSizeCurrent, nd, Chain);
                else
                    sampleSizeCurrent = sum(ChainWeight);
                    [CovMatUpperCurrent, MeanCurrent] = getWeiSamCovUppMeanTrans(chainSize, sampleSizeCurrent, nd, Chain, ChainWeight);
                end

                %***********************************************************************************************************************

                % combine old and new covariance matrices if both exist

                if self.Global.mv_sampleSizeOld_save == 0   % blockMergeCovMat

                    % There is no prior old Covariance matrix to combine with the new one from the new chain

                    self.Global.mv_MeanOld_save(1 : nd) = MeanCurrent;
                    self.Global.mv_sampleSizeOld_save   = sampleSizeCurrent;

                    % copy and then scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor

                    CovMatUpperOld      = self.Global.comv_covMat(:,:,1);
                    self.Global.comv_covMat(:,:,1)  = CovMatUpperCurrent * self.Global.mc_defaultScaleFactorSq;

                else % blockMergeCovMat

                    % first scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor

                    CovMatUpperOld      = self.Global.comv_covMat(:,:,1);
                    CovMatUpperCurrent  = CovMatUpperCurrent * self.Global.mc_defaultScaleFactorSq;

                    % now combine it with the old covariance matrix

                    [MeanNew, self.Global.comv_covMat(:,:,1)] = combineMeanCovUpper ( nd                                    ...
                                                                        , self.Global.mv_sampleSizeOld_save                 ...
                                                                        , self.Global.mv_MeanOld_save                       ...
                                                                        , CovMatUpperOld                        ...
                                                                        , sampleSizeCurrent                     ...
                                                                        , MeanCurrent                           ...
                                                                        , CovMatUpperCurrent                    ...
                                                                        ) ;

                    self.Global.mv_MeanOld_save(1:nd) = MeanNew;

                end % blockMergeCovMat

                %***********************************************************************************************************************

                % now get the Cholesky factorization

                self.Global.comv_chol = chol(self.Global.comv_covMat(:, :, 1), 'lower');

                % if comv_chol(1, 1, 1) > 0  % blockPosDefCheck
                
                    ... singularityOccurred         = false;
                    samplerUpdateSucceeded      = true;
                    hellingerComputationNeeded  = true;
                    self.Global.mv_sampleSizeOld_save       = self.Global.mv_sampleSizeOld_save + sampleSizeCurrent;
                
                % else % blockPosDefCheck
                % 
                %     % it may be a good idea to add a warning message printed out here for the singularity occurrence
                % 
                %     singularityOccurred         = true;
                %     samplerUpdateSucceeded      = false;
                %     hellingerComputationNeeded  = false;
                %     mv_sampleSizeOld_save       = sampleSizeOld;
                % 
                %     % recover the old upper covariance matrix
                % 
                %     comv_covMat(:, :, 1) = CovMatUpperOld;
                % 
                %   % for j = 2 : nd+1
                %       % for i = 1 : j-1
                %           % comv_CholDiagLower(i,j,1) = CovMatUpperOld(i,j);
                %       % end
                %   % end
                % 
                %     % ensure the old Cholesky factorization can be recovered
                % 
                %     comv_chol = chol(comv_covMat(:,:,1), 'lower');
                %   % getCholeskyFactor(nd, comv_CholDiagLower(1:nd, 2:nd+1, 1), comv_CholDiagLower(1:nd,1,1));
                %     if comv_chol(1,1,1) < 0
                %         mv_Err.msg  = FUNCTION_NAME                                                                                                     ...
                %                     + ": Error occurred while attempting to compute the Cholesky factorization of the "                                 ...
                %                     + "covariance matrix of the proposal distribution of " + mc_methodName + "'s sampler. "                             ...
                %                     + "This is highly unusual, and can be indicative of some major underlying problems. "                               ...
                %                     + "It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations." + newline...
                %                     + "In such case, consider increasing the value of the input variable adaptiveUpdatePeriod." + newline               ...
                %                     + "Restarting the simulation might resolve the error."                                                              ...
                %                     ;
                %         mv_Err.abort(mc_methodBrand, newline, mc_logFileUnit);
                %     end
                % end % blockPosDefCheck

                if self.Global.mc_scalingRequested, scalingNeeded = true; end

            else % blockSufficientSampleSizeCheck

                % singularity has occurred. If the first covariance merging has not occurred yet, set the scaling factor appropriately to shrink the covariance matrix.

                samplerUpdateSucceeded = false;
                if (self.Global.mv_sampleSizeOld_save == 0) || self.Global.mc_scalingRequested
                    scalingNeeded               = true;
                    hellingerComputationNeeded  = true;
                    % save the old covariance matrix for the computation of adaptation measure
                    CovMatUpperOld = self.Global.comv_covMat;
                else
                    hellingerComputationNeeded  = false;
                end

            end % blockSufficientSampleSizeCheck

            % adjust the scale of the covariance matrix and the Cholesky factor, if needed

            if class(meanAccRateSinceStart) == "char", meanAccRateSinceStart = str2num(meanAccRateSinceStart); end

            if scalingNeeded

                self.Global.mv_adaptiveScaleFactorSq_save   = (meanAccRateSinceStart / self.Global.mc_targetAcceptanceRate)^self.Global.mc_ndimInverse;
                adaptiveScaleFactor             = sqrt(self.Global.mv_adaptiveScaleFactorSq_save);

                % update the Cholesky diagonal elements
                self.Global.comv_chol(:,:,1)    = self.Global.comv_chol(:,:,1) * adaptiveScaleFactor;

                % update covariance matrix
                self.Global.comv_covMat         = self.Global.comv_covMat * self.Global.mv_adaptiveScaleFactorSq_save;

            end

            % compute the adaptivity only if any updates has occurred

            if hellingerComputationNeeded
                logSqrtDetNew = sum(log(diag(self.Global.comv_chol(:,:,1))));
                CovMatUpperCurrent  = 0.5 * (self.Global.comv_covMat(:,:,1) + CovMatUpperOld);

                try
                    logSqrtDetSum   = 0.5 * log(det(CovMatUpperCurrent(:,:,1)));
                catch exception
                    self.Global.mv_Err.msg  = FUNCTION_NAME                                                                                         ...
                                + ": Error occurred while computing the Cholesky factorization of "                                                 ...
                                + "a matrix needed for the computation of the Adaptation measure. "                                                 ...
                                + "Such error is highly unusual, and requires an in depth investigation of the case. "                              ...
                                + "It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations." + newline...
                                + "In such case, consider increasing the value of the input variable adaptiveUpdatePeriod." + newline               ...
                                + "Restarting the simulation might resolve the error."                                                              ...
                                ;
                    return
                end

              % getLogSqrtDetPosDefMat(nd, CovMatUpperCurrent, logSqrtDetSum, singularityOccurred);
              % if singularityOccurred
                  % mv_Err.msg  = FUNCTION_NAME                                                                                                     ...
                              % + ": Error occurred while computing the Cholesky factorization of "                                                 ...
                              % + "a matrix needed for the computation of the Adaptation measure. "                                                 ...
                              % + "Such error is highly unusual, and requires an in depth investigation of the case. "                              ...
                              % + "It may also be due to a runtime computational glitch, in particular, for high-dimensional simulations." + newline...
                              % + "In such case, consider increasing the value of the input variable adaptiveUpdatePeriod." + newline               ...
                              % + "Restarting the simulation might resolve the error."                                                              ...
                              % ;
                  % mv_Err.abort(mc_methodBrand, newline, mc_logFileUnit);
                  % return
              % end
                adaptationMeasure = sqrt(1.0 - exp(0.5 * (self.Global.mv_logSqrtDetOld_save + logSqrtDetNew) - logSqrtDetSum));
                if adaptationMeasure < 0.0
                    self.Global.mv_Err.msg  = self.Global.mc_negativeAdaptationMeasureMsg + num2str(adaptationMeasure);
                    self.Global.mv_Err.warn();
                    adaptationMeasure         = 0.0;
                end
                % update the higher-stage delayed-rejection Cholesky Lower matrices
                if self.Global.mc_delayedRejectionRequested, self.updateDelRejCholDiagLower(); end
            end

        end % doAdaptation

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function updateDelRejCholDiagLower(self)

%            global  mc_DelayedRejectionCount            ...
%                    comv_chol                           ...
%                    mc_DelayedRejectionScaleFactorVec

            % update the Cholesky factor of the delayed-rejection-stage proposal distributions
            for istage = 2 : self.Global.mc_DelayedRejectionCount + 1
                self.Global.comv_chol(:,:,istage) = self.Global.comv_chol(:,:,istage-1) * self.Global.mc_DelayedRejectionScaleFactorVec(istage-1);
            end
            % There is no need to check for positive-definiteness of the comv_CholDiagLower, it's already checked on the first image.

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function writeRestartFile(self)
            
            k = 1;
            for j = 1 : self.Global.mc_ndim
                for i = 1 : j
                    self.Global.lower_comv_covMat(k) = self.Global.comv_covMat(i,j);
                    k = k + 1;
                end
            end

            fprintf(self.Global.mc_restartFileUnit  , self.Global.mc_restartFileFormat          ...
                                                    , self.Global.mv_sampleSizeOld_save         ...
                                                    , self.Global.mv_logSqrtDetOld_save         ...
                                                    , self.Global.mv_adaptiveScaleFactorSq_save ...
                                                    , self.Global.mv_MeanOld_save.'             ...
                                                    ..., self.Global.comv_chol(:,:,1).'         ... not needed anymore since lower comv_covMat is sufficient
                                                    ..., self.Global.comv_covMat(:,:,1).'       ...
                                                    , self.Global.lower_comv_covMat(1,:)        ...
                                                    );

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function readRestartFile(self)

            for i = 1 : 9 + self.Global.mc_ndim + 2*(self.Global.mc_ndim)^2
                fgets(self.Global.mc_restartFileUnit);
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end