classdef ParaDRAMProposalSymmetric_class < handle

    properties (Constant)
        SUB_CLASS_NAME      = "@ParaDRAMProposalSymmetric_class"
    end

    properties
       %AccRate             = struct('sumUpToLastUpdate', [], 'target', [])
        Err                 = Err_class()
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

            global  comv_chol                           ...
                    comv_covMat                         ...
                    mc_Image                            ...
                    mc_ndim                             ...
                    mc_logFileUnit                      ...
                    mc_restartFileUnit                  ...
                    mc_scalingRequested                 ...
                    mc_defaultScaleFactorSq             ...
                    mc_DelayedRejectionCount            ...
                    mc_MaxNumDomainCheckToWarn          ...
                    mc_MaxNumDomainCheckToStop          ...
                    mc_delayedRejectionRequested        ...
                    mc_ndimInverse                      ...
                    mc_targetAcceptanceRate             ...
                    mc_DelayedRejectionScaleFactorVec   ...
                    mc_DomainLowerLimitVec              ...
                    mc_DomainUpperLimitVec              ...
                    mc_MaxNumDomainCheckToWarnMsg       ...
                    mc_MaxNumDomainCheckToStopMsg       ...
                    mc_negativeHellingerDistSqMsg       ...
                    mc_restartFileFormat                ...
                    mc_methodBrand                      ...
                    mc_methodName                       ...
                    mc_isNormal                         ...
                    ...the following had to be defined globally for the sake of restart file generation
                    mv_MeanOld_save                     ...
                    mv_logSqrtDetOld_save               ...
                    mv_adaptiveScaleFactorSq_save       ...
                    mv_sampleSizeOld_save               ...
                    mv_Err

            %***********************************************************************************************************************
            % setup sampler update global save variables
            %***********************************************************************************************************************

            mv_MeanOld_save(1:ndim)         = Constants.NULL_RK;
            mv_logSqrtDetOld_save           = Constants.NULL_RK;
            mv_sampleSizeOld_save           = 0;
            mv_adaptiveScaleFactorSq_save   = 1;
            mv_Err                          = Err_class();

            %***********************************************************************************************************************
            % setup general ParaDRAMProposalSymmetric specifications
            %***********************************************************************************************************************

            mc_ndim                             = ndim;
            mc_DomainLowerLimitVec              = PD.SpecBase.domainLowerLimitVec.Val;
            mc_DomainUpperLimitVec              = PD.SpecBase.domainUpperLimitVec.Val;
            mc_DelayedRejectionScaleFactorVec   = PD.SpecDRAM.delayedRejectionScaleFactorVec.Val;
            mc_isNormal                         = PD.SpecDRAM.proposalModel.isNormal;
            mc_Image                            = PD.Image;
            mc_methodName                       = PD.name;
            mc_methodBrand                      = PD.brand;
            mc_logFileUnit                      = PD.LogFile.unit;
            mc_restartFileUnit                  = PD.RestartFile.unit;
            mc_restartFileFormat                = PD.RestartFile.format;
            mc_defaultScaleFactorSq             = PD.SpecDRAM.scaleFactor.val^2;
           %self.AccRate.sumUpToLastUpdate      = 0;
            mc_maxNumDomainCheckToWarn          = PD.SpecBase.maxNumDomainCheckToWarn.val;
            mc_maxNumDomainCheckToStop          = PD.SpecBase.maxNumDomainCheckToStop.val;
            mc_DelayedRejectionCount            = PD.SpecDRAM.delayedRejectionCount.val;
            mc_delayedRejectionRequested        = mc_DelayedRejectionCount > 0;
            mc_scalingRequested                 = PD.SpecBase.targetAcceptanceRate.scalingRequested;

            if mc_scalingRequested
                mc_targetAcceptanceRate         = PD.SpecBase.targetAcceptanceRate.val;
            else
                mc_targetAcceptanceRate         = 1.0;    % 0.234
            end

            mc_ndimInverse                      = 1.0 / ndim;

            %***********************************************************************************************************************
            % setup ProposalSymmetric specifications
            %***********************************************************************************************************************

            % setup covariance matrix

            comv_covMat = zeros(ndim, ndim, mc_DelayedRejectionCount + 1); % xxx: define new variable 
            comv_chol   = zeros(ndim, ndim, mc_DelayedRejectionCount + 1);

            if PD.SpecDRAM.proposalStartCovMat.isPresent
                comv_covMat(:, :, 1)  = PD.SpecDRAM.proposalStartCovMat.Val;
            else
                comv_covMat(:, :, 1)  = getCovMatFromCorMat ( ndim                                  ...
                                                            , PD.SpecDRAM.proposalStartStdVec.Val   ...
                                                            , PD.SpecDRAM.proposalStartCorMat.Val   ...
                                                            ) ;
            end

            % Now scale the covariance matrix

            comv_covMat(:,:,1) = comv_covMat(:,:,1) * mc_defaultScaleFactorSq;

            % Now get the Cholesky Factor of the Covariance Matrix. Lower comv_CholDiagLower will be the CholFac

            comv_chol(:,:,1)    = chol(comv_covMat(:, :, 1), 'lower'); % xxx check what happens if chol fails

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

            if mc_delayedRejectionRequested, ParaDRAMProposalSymmetric_class.updateDelRejCholDiagLower(); end

            % This will be used for Domain boundary checking during the simulation

            mc_negativeHellingerDistSqMsg   = self.SUB_CLASS_NAME + "@doAdaptation(): Non-positive adaptation measure detected, "           ...
                                            + "possibly due to round-off error: ";

            mc_MaxNumDomainCheckToWarnMsg   = self.SUB_CLASS_NAME + "@getNew(): " + num2str(mc_MaxNumDomainCheckToWarn) + " proposals were "...
                                            + "drawn out of the objective function's Domain without any acceptance.";

            mc_MaxNumDomainCheckToStopMsg   = self.SUB_CLASS_NAME + "@getNew(): " + num2str(mc_MaxNumDomainCheckToStop)                     ...
                                            + " proposals were drawn out of the objective function's Domain. As per the value set "         ...
                                            + "for the simulation specification variable 'maxNumDomainCheckToStop', " + mc_methodName       ...
                                            + " will abort now.";

        end % ParaDRAMProposalSymmetric_class

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        function StateNew = getNew  ( nd            ...
                                    , counterDRS    ...
                                    , StateOld      ...
                                    )
            StateOld = transpose(StateOld);
    
            FUNCTION_NAME = ParaDRAMProposalSymmetric_class.SUB_CLASS_NAME + "@getNew()";

            global  comv_chol                       ...
                    comv_covMat                     ...
                    mc_logFileUnit                  ...
                    mc_MaxNumDomainCheckToWarn      ...
                    mc_MaxNumDomainCheckToStop      ...
                    mc_DomainLowerLimitVec          ...
                    mc_DomainUpperLimitVec          ...
                    mc_MaxNumDomainCheckToWarnMsg   ...
                    mc_MaxNumDomainCheckToStopMsg   ...
                    mc_methodBrand                  ...
                    mc_isNormal                     ...
                    mv_Err
                    
            mv_Err.prefix       = mc_methodBrand;
            mv_Err.outputUnit   = mc_logFileUnit;

            domainCheckCounter = 0;

            if mc_isNormal

                while true  % loopBoundaryCheckNormal: Check for the support Region consistency:

%                    StateNew = mvnrnd(StateOld, comv_chol(:,:,counterDRS+1));

                    StateNew = getRandMVN   ( nd                            ...
                                            , StateOld                      ...
                                            , comv_chol (:,:,counterDRS+1)  ...
                                            ) ;

                    if any(StateNew <= mc_DomainLowerLimitVec) || any(StateNew >= mc_DomainUpperLimitVec)
                        domainCheckCounter = domainCheckCounter + 1;
                        if domainCheckCounter == mc_MaxNumDomainCheckToWarn
                            mv_Err.msg      = mc_MaxNumDomainCheckToWarnMsg;
                            mv_Err.warn();
                        end
                        if domainCheckCounter == mc_MaxNumDomainCheckToStop
                            mv_Err.msg      = mc_MaxNumDomainCheckToStopMsg;
                            mv_Err.abort();
                        end
                        continue
                    end
                    break
                end % loopBoundaryCheckNormal

            else % if (mc_isUniform)

                while true  % loopBoundaryCheckUniform

                    StateNew = getRandMVU   ( nd                            ...
                                            , StateOld                      ...
                                            , comv_chol (:,:,counterDRS+1)  ...
                                            ) ;

                    if any(StateNew <= mc_DomainLowerLimitVec) || any(StateNew >= mc_DomainUpperLimitVec)
                        domainCheckCounter  = domainCheckCounter + 1;
                        if domainCheckCounter == mc_MaxNumDomainCheckToWarn
                            mv_Err.msg      = mc_MaxNumDomainCheckToWarnMsg;
                            mv_Err.warn();
                        end
                        if domainCheckCounter == mc_MaxNumDomainCheckToStop
                            mv_Err.msg      = mc_MaxNumDomainCheckToStopMsg;
                            mv_Err.abort();
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
        function hellingerDistSq =  doAutoTune  ( hellingerDistSq   ...
                                                , AutoTuneScaleSq   ...
                                                )

            FUNCTION_NAME           = SUB_CLASS_NAME + "@doAutoTune()";
            
            mv_Err.prefix           = mc_methodBrand;
            mv_Err.outputUnit       = mc_logFileUnit;
            
            CovMatUpperOld          = comv_CholDiagLower(1:mc_ndim, 2:mc_ndim+1, 1);
            mv_logSqrtDetOld_save   = sum(log(comv_CholDiagLower(1:mc_ndim, 1, 1)));

            if AutoTuneScaleSq(1) == 0
                comv_CholDiagLower(1,2,1) = 0.25 * comv_CholDiagLower(1,2,1);
                comv_CholDiagLower(1,1,1) = sqrt(comv_CholDiagLower(1,2,1));
            else
                comv_CholDiagLower(1,2,1) = AutoTuneScaleSq(1);
                comv_CholDiagLower(1,1,1) = sqrt(AutoTuneScaleSq(1));
            end

            % compute the hellingerDistSq (adaptivity)
            logSqrtDetNew = sum(log(comv_CholDiagLower(1:mc_ndim, 1, 1)));
            CovMatUpperCurrent = 0.5 * (comv_CholDiagLower(1:mc_ndim, 2:mc_ndim+1, 1) + CovMatUpperOld);
            [logSqrtDetSum, singularityOccurred] = getLogSqrtDetPosDefMat(1, CovMatUpperCurrent);
            if singularityOccurred
                mv_Err.msg  = FUNCTION_NAME                                                                             ...
                            + ": Error occurred while computing the Cholesky factorization of "                         ...
                            + "a matrix needed for the computation of the proposal distribution's adaptation measure. " ...
                            + "Such error is highly unusual, and requires an in depth investigation of the case. "      ...
                            + "Restarting the simulation might resolve the error."                                      ...
                            ;
                mv_Err.abort();
                return
            end
            hellingerDistSq = 1 - exp(0.5 * (mv_logSqrtDetOld_save + logSqrtDetNew) - logSqrtDetSum);

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [samplerUpdateSucceeded, hellingerDistSq]  = doAdaptation  ( nd                        ...
                                                                            , chainSize                 ...
                                                                            , Chain                     ...
                                                                            , ChainWeight               ...
                                                                            , samplerUpdateIsGreedy     ...
                                                                            , meanAccRateSinceStart     ...
                                                                            , hellingerDistSq           ...
                                                                            )
            FUNCTION_NAME = ParaDRAMProposalSymmetric_class.SUB_CLASS_NAME + "@doAdaptation()";

            global  comv_chol                           ...
                    comv_covMat                         ...
                    mc_Image                            ...
                    mc_ndim                             ...
                    mc_logFileUnit                      ...
                    mc_restartFileUnit                  ...
                    mc_scalingRequested                 ...
                    mc_defaultScaleFactorSq             ...
                    mc_DelayedRejectionCount            ...
                    mc_MaxNumDomainCheckToWarn          ...
                    mc_MaxNumDomainCheckToStop          ...
                    mc_delayedRejectionRequested        ...
                    mc_ndimInverse                      ...
                    mc_targetAcceptanceRate             ...
                    mc_DelayedRejectionScaleFactorVec   ...
                    mc_DomainLowerLimitVec              ...
                    mc_DomainUpperLimitVec              ...
                    mc_MaxNumDomainCheckToWarnMsg       ...
                    mc_MaxNumDomainCheckToStopMsg       ...
                    mc_negativeHellingerDistSqMsg       ...
                    mc_restartFileFormat                ...
                    mc_methodBrand                      ...
                    mc_methodName                       ...
                    mc_isNormal                         ...
                    ...the following had to be defined globally for the sake of restart file generation
                    mv_MeanOld_save                     ...
                    mv_logSqrtDetOld_save               ...
                    mv_adaptiveScaleFactorSq_save       ...
                    mv_sampleSizeOld_save               ...
                    mv_Err

            scalingNeeded           = false;
            sampleSizeOld           = mv_sampleSizeOld_save;  % this is kept only for restoration of mv_sampleSizeOld_save, if needed.
            mv_logSqrtDetOld_save   = sum(log(diag(comv_chol(:,:,1))));

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

                if mv_sampleSizeOld_save == 0   % blockMergeCovMat

                    % There is no prior old Covariance matrix to combine with the new one from the new chain

                    mv_MeanOld_save(1 : nd) = MeanCurrent;
                    mv_sampleSizeOld_save   = sampleSizeCurrent;

                    % copy and then scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor

                    CovMatUpperOld      = comv_covMat(:,:,1);
                    comv_covMat(:,:,1)  = CovMatUpperCurrent * mc_defaultScaleFactorSq;

                else % blockMergeCovMat

                    % first scale the new covariance matrix by the default scale factor, which will be then used to get the Cholesky Factor

                    CovMatUpperOld      = comv_covMat(:,:,1);
                    CovMatUpperCurrent  = CovMatUpperCurrent * mc_defaultScaleFactorSq;

                    % now combine it with the old covariance matrix

                    [MeanNew, comv_covMat(:,:,1)] = combineMeanCovUpper ( nd                                    ...
                                                                        , mv_sampleSizeOld_save                 ...
                                                                        , mv_MeanOld_save                       ...
                                                                        , CovMatUpperOld                        ...
                                                                        , sampleSizeCurrent                     ...
                                                                        , MeanCurrent                           ...
                                                                        , CovMatUpperCurrent                    ...
                                                                        ) ;

                    mv_MeanOld_save(1:nd) = MeanNew;

                end % blockMergeCovMat

                %***********************************************************************************************************************

                % now get the Cholesky factorization

                comv_chol = chol(comv_covMat(:, :, 1), 'lower');

                % if comv_chol(1, 1, 1) > 0  % blockPosDefCheck
                
                    singularityOccurred         = false;
                    samplerUpdateSucceeded      = true;
                    hellingerComputationNeeded  = true;
                    mv_sampleSizeOld_save       = mv_sampleSizeOld_save + sampleSizeCurrent;
                
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

                if mc_scalingRequested, scalingNeeded = true; end

            else % blockSufficientSampleSizeCheck

                % singularity has occurred. If the first covariance merging has not occurred yet, set the scaling factor appropriately to shrink the covariance matrix.

                samplerUpdateSucceeded = false;
                if (mv_sampleSizeOld_save == 0) || mc_scalingRequested
                    scalingNeeded               = true;
                    hellingerComputationNeeded  = true;
                    % save the old covariance matrix for the computation of adaptation measure
                    CovMatUpperOld = comv_covMat;
                else
                    hellingerComputationNeeded  = false;
                end

            end % blockSufficientSampleSizeCheck

            % adjust the scale of the covariance matrix and the Cholesky factor, if needed

            if class(meanAccRateSinceStart) == "char", meanAccRateSinceStart = str2num(meanAccRateSinceStart); end
            if scalingNeeded

                mv_adaptiveScaleFactorSq_save   = (meanAccRateSinceStart / mc_targetAcceptanceRate)^mc_ndimInverse;
                adaptiveScaleFactor             = sqrt(mv_adaptiveScaleFactorSq_save);
              
                % update the Cholesky diagonal elements
                comv_chol(:,:,1)    = comv_chol(:,:,1) * adaptiveScaleFactor;
                
                % update covariance matrix
                comv_covMat         = comv_covMat * mv_adaptiveScaleFactorSq_save;
              
            end

            % compute the adaptivity only if any updates has occurred

            if hellingerComputationNeeded
                logSqrtDetNew = sum(log(diag(comv_chol)));
                CovMatUpperCurrent  = 0.5 * (comv_covMat(:,:,1) + CovMatUpperOld);
              
                try
                    logSqrtDetSum   = 0.5 * log(det(CovMatUpperCurrent(:,:,1)));
                catch exception
                    mv_Err.msg  = FUNCTION_NAME                                                                                                     ...
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
                hellingerDistSq = 1.0 - exp(0.5 * (mv_logSqrtDetOld_save + logSqrtDetNew) - logSqrtDetSum);
                if hellingerDistSq < 0.0
                    Err_class.warn(mc_negativeHellingerDistSqMsg + num2str(hellingerDistSq), mc_methodBrand, [], mc_logFileUnit);
                    hellingerDistSq = 0.0;
                end
                % update the higher-stage delayed-rejection Cholesky Lower matrices
                if mc_delayedRejectionRequested, ParaDRAMProposalSymmetric_class.updateDelRejCholDiagLower(); end
            end

        end % doAdaptation

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function updateDelRejCholDiagLower()

            global  mc_DelayedRejectionCount            ...
                    comv_chol                           ...
                    mc_DelayedRejectionScaleFactorVec

            % update the Cholesky factor of the delayed-rejection-stage proposal distributions
            for istage = 2 : mc_DelayedRejectionCount + 1

                comv_chol(:,:,istage) = comv_chol(:,:,istage-1) * mc_DelayedRejectionScaleFactorVec(istage-1);

            end
            % There is no need to check for positive-definiteness of the comv_CholDiagLower, it's already checked on the first image.
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function writeRestartFile()

            global  comv_chol                       ...
                    comv_covMat                     ...
                    mc_restartFileUnit              ...
                    mc_restartFileFormat            ...
                    mv_sampleSizeOld_save           ...
                    mv_logSqrtDetOld_save           ...
                    mv_adaptiveScaleFactorSq_save   ...
                    mv_MeanOld_save                 ...
                    mc_ndim

            fprintf(mc_restartFileUnit, "\n%s"  , "sampleSizeOld"                           ...
                                                , num2str(mv_sampleSizeOld_save)            ...
                                                , "logSqrtDetOld"                           ...
                                                , num2str(mv_logSqrtDetOld_save)            ...
                                                , "adaptiveScaleFactorSq"                   ...
                                                , num2str(mv_adaptiveScaleFactorSq_save)    ...
                                                ) ;

            fprintf(mc_restartFileUnit, "\n%s", "MeanOld");
            fprintf(mc_restartFileUnit, "\n%d", mv_MeanOld_save.');

            fprintf(mc_restartFileUnit, "\n%s", "comv_chol");
            fprintf(mc_restartFileUnit, "\n%d", comv_chol(:,:,1).');
            
            fprintf(mc_restartFileUnit, "\n%s", "comv_covMat");
            fprintf(mc_restartFileUnit, "\n%d", comv_covMat(:,:,1).');
            
            fprintf(mc_restartFileUnit, "\n");

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function readRestartFile()

            global mc_restartFileUnit mc_ndim

            for i = 1 : 8 + mc_ndim * (mc_ndim + 2)
                fgets(mc_restartFileUnit);
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end