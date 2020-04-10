classdef SpecDRAM_class < handle

    properties (Constant)
    end

    properties
        scaleFactor                     = []
        proposalModel                   = []
        proposalStartCovMat             = []
        proposalStartCorMat             = []
        proposalStartStdVec             = []
        adaptiveUpdateCount             = []
        adaptiveUpdatePeriod            = []
        greedyAdaptationCount           = []
        delayedRejectionCount           = []
        burninAdaptationMeasure         = []
        delayedRejectionScaleFactorVec  = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_class( nd, methodName)
            self.scaleFactor                    = SpecDRAM_ScaleFactor_class                    (nd,methodName);
            self.proposalModel                  = SpecDRAM_ProposalModel_class                  ( );
            self.proposalStartCovMat            = SpecDRAM_ProposalStartCovMat_class            (nd,methodName);
            self.proposalStartCorMat            = SpecDRAM_ProposalStartCorMat_class            (nd,methodName);
            self.proposalStartStdVec            = SpecDRAM_ProposalStartStdVec_class            (nd,methodName);
            self.adaptiveUpdatePeriod           = SpecDRAM_AdaptiveUpdatePeriod_class           (nd,methodName);
            % ATTN: AdaptiveUpdateCount has to be constructed after AdaptiveUpdatePeriod. It depends on it.
            self.adaptiveUpdateCount            = SpecDRAM_AdaptiveUpdateCount_class            (methodName); % ,chainSizeDef, self.AdaptiveUpdatePeriod.def)
            self.delayedRejectionCount          = SpecDRAM_DelayedRejectionCount_class          (methodName);
            self.burninAdaptationMeasure        = SpecDRAM_BurninAdaptationMeasure_class        (methodName);
            self.greedyAdaptationCount          = SpecDRAM_GreedyAdaptationCount_class          (methodName);
            self.delayedRejectionScaleFactorVec = SpecDRAM_DelayedRejectionScaleFactorVec_class (nd,methodName);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function setFromInputArgs   ( self                              ...
                                    ... % input arguments to the specific ParaDRAM routine
                                    , scaleFactor                       ...
                                    , proposalModel                     ...
                                    , proposalStartCovMat               ...
                                    , proposalStartCorMat               ...
                                    , proposalStartStdVec               ...
                                    , adaptiveUpdateCount               ...
                                    , adaptiveUpdatePeriod              ...
                                    , greedyAdaptationCount             ...
                                    , delayedRejectionCount             ...
                                    , burninAdaptationMeasure           ...
                                    , delayedRejectionScaleFactorVec    ...
                                    )

            self.scaleFactor.                       set(scaleFactor);
            self.proposalModel.                     set(proposalModel);
            self.proposalStartCovMat.               set(proposalStartCovMat);
            self.proposalStartCorMat.               set(proposalStartCorMat);
            self.proposalStartStdVec.               set(proposalStartStdVec);
            self.adaptiveUpdateCount.               set(adaptiveUpdateCount);
            self.adaptiveUpdatePeriod.              set(adaptiveUpdatePeriod);
            self.greedyAdaptationCount.             set(greedyAdaptationCount);
            self.delayedRejectionCount.             set(delayedRejectionCount);
            self.burninAdaptationMeasure.           set(burninAdaptationMeasure);
            self.delayedRejectionScaleFactorVec.    set(delayedRejectionScaleFactorVec, self.delayedRejectionCount.val);

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function reportValues(self, prefix, outputUnit, methodName, splashModeRequested)

            formatVal = Decoration_class.TAB + Decoration_class.TAB;

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "adaptiveUpdatePeriod" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.adaptiveUpdatePeriod.val) + "\n");
            if splashModeRequested, Err_class.note(self.adaptiveUpdatePeriod.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "adaptiveUpdateCount" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.adaptiveUpdateCount.val) + "\n");
            if splashModeRequested, Err_class.note(self.adaptiveUpdateCount.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "greedyAdaptationCount" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.greedyAdaptationCount.val) + "\n");
            if splashModeRequested, Err_class.note(self.greedyAdaptationCount.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "burninAdaptationMeasure" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.burninAdaptationMeasure.val) + "\n");
            if splashModeRequested, Err_class.note(self.burninAdaptationMeasure.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "delayedRejectionCount" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.delayedRejectionCount.val) + "\n");
            if splashModeRequested, Err_class.note(self.delayedRejectionCount.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "delayedRejectionScaleFactorVec" + "\n\n");
            if isempty(self.delayedRejectionScaleFactorVec.Val)
                fprintf(outputUnit, formatVal + "UNDEFINED" + "\n");
            else
                for i = 1 : length(self.delayedRejectionScaleFactorVec.Val)
                    fprintf(outputUnit, formatVal + num2str(self.delayedRejectionScaleFactorVec.Val(i)) + "\n");
                end
            end
            if splashModeRequested, Err_class.note(self.delayedRejectionScaleFactorVec.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "scaleFactor" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.scaleFactor.str) + "\n");
            if splashModeRequested, Err_class.note(self.scaleFactor.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            %***********************************************************************************************************************
            % proposal distribution
            %***********************************************************************************************************************

            Decoration = Decoration_class([],[],[],[]);

            Decoration.writeDecoratedText(" " + newline + methodName + " proposal specifications" + newline, [], [], [], [], 1, 1, outputUnit, "\n");


            fprintf(outputUnit,  "\n" + "proposalModel" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.proposalModel.val) + "\n");
            if splashModeRequested, Err_class.note(self.proposalModel.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            ndim = length(self.proposalStartCovMat.Val(:,1));

            fprintf(outputUnit,  "\n" + "proposalStartCovMat" + "\n\n");
            if self.proposalStartCovMat.isPresent
                % User has provided the Start Covariance Matrix
                for i = 1 : ndim
                    Row = self.proposalStartCovMat.Val(i,:);
                    fprintf(outputUnit, formatVal + num2str(Row) + "\n");
                end
            else
                % User has not provided the Start Covariance Matrix
                Err_class.informUser("UNDEFINED. It will be constructed from the Correlation Matrix (ProposalStartCorMat)"    ...
                                        + "and the Standard Deviation vector (ProposalStartStdVec)."                          ...
                                        , Decoration_class.TAB + Decoration_class.TAB, "\n", outputUnit, [], [], [], 0);
            end
            if splashModeRequested, Err_class.note(self.proposalStartCovMat.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "proposalStartCorMat" + "\n\n");
            for i = 1 : ndim
                Row = self.proposalStartCorMat.Val(i,:);
                fprintf(outputUnit, formatVal + num2str(Row) + "\n");
            end
            if splashModeRequested, Err_class.note(self.proposalStartCorMat.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "proposalStartStdVec" + "\n\n");
            for i = 1 : ndim
                fprintf(outputUnit, formatVal + num2str(self.proposalStartStdVec.Val(i)) + "\n");
            end
            if splashModeRequested, Err_class.note(self.proposalStartStdVec.desc, prefix, "\n", outputUnit, [], []); end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, nd)
            self.scaleFactor.                       checkForSanity (Err, methodName)                                    ;
            self.proposalModel.                     checkForSanity (Err, methodName)                                    ;
            self.adaptiveUpdateCount.               checkForSanity (Err, methodName)                                    ;
            self.adaptiveUpdatePeriod.              checkForSanity (Err, methodName)                                    ;
            self.greedyAdaptationCount.             checkForSanity (Err, methodName)                                    ;
            self.delayedRejectionCount.             checkForSanity (Err, methodName)                                    ;
            self.burninAdaptationMeasure.           checkForSanity (Err, methodName)                                    ;
            self.proposalStartCovMat.               checkForSanity (Err, methodName, nd)                                ;
            self.proposalStartCorMat.               checkForSanity (Err, methodName, nd)                                ;
            self.proposalStartStdVec.               checkForSanity (Err, methodName, nd)                                ;
            self.delayedRejectionScaleFactorVec.    checkForSanity (Err, methodName, self.delayedRejectionCount.val)    ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end