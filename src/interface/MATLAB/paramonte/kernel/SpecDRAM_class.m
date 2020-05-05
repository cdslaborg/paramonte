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

            formatVal           = Decoration_class.TAB + Decoration_class.TAB;
            
            Err                 = Err_class();
            Err.prefix          = prefix;
            Err.outputUnit      = outputUnit;
            Err.resetEnabled    = false;

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "adaptiveUpdatePeriod" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.adaptiveUpdatePeriod.val) + "\n");
            Err.msg             = self.adaptiveUpdatePeriod.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "adaptiveUpdateCount" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.adaptiveUpdateCount.val) + "\n");
            Err.msg             = self.adaptiveUpdateCount.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "greedyAdaptationCount" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.greedyAdaptationCount.val) + "\n");
            Err.msg             = self.greedyAdaptationCount.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "burninAdaptationMeasure" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.burninAdaptationMeasure.val) + "\n");
            Err.msg             = self.burninAdaptationMeasure.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "delayedRejectionCount" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.delayedRejectionCount.val) + "\n");
            Err.msg             = self.delayedRejectionCount.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "delayedRejectionScaleFactorVec" + "\n\n");
            if isempty(self.delayedRejectionScaleFactorVec.Val)
                fprintf(outputUnit, formatVal + "UNDEFINED" + "\n");
            else
                for i = 1 : length(self.delayedRejectionScaleFactorVec.Val)
                    fprintf(outputUnit, formatVal + num2str(self.delayedRejectionScaleFactorVec.Val(i)) + "\n");
                end
            end
            Err.msg             = self.delayedRejectionScaleFactorVec.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "scaleFactor" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.scaleFactor.str) + "\n");
            Err.msg             = self.scaleFactor.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            %***********************************************************************************************************************
            % proposal distribution
            %***********************************************************************************************************************

            Decoration          = Decoration_class([],[],[],[]);

            Decoration.writeDecoratedText(" " + newline + methodName + " proposal specifications" + newline, [], [], [], [], 1, 1, outputUnit, "\n");


            fprintf(outputUnit,  "\n" + "proposalModel" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.proposalModel.val) + "\n");
            Err.msg             = self.proposalModel.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            ndim                = length(self.proposalStartCovMat.Val(:,1));

            fprintf(outputUnit,  "\n" + "proposalStartCovMat" + "\n\n");
            if self.proposalStartCovMat.isPresent
                % User has provided the Start Covariance Matrix
                for i = 1 : ndim
                    Row = self.proposalStartCovMat.Val(i,:);
                    fprintf(outputUnit, formatVal + num2str(Row) + "\n");
                end
            else
                % User has not provided the Start Covariance Matrix
                Err.fullprefix  = formatVal;
                Err.msg         = "UNDEFINED. It will be constructed from the Correlation Matrix (ProposalStartCorMat)"        ...
                                + "and the Standard Deviation vector (ProposalStartStdVec).";
                Err.informUser();
            end
            Err.msg             = self.proposalStartCovMat.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "proposalStartCorMat" + "\n\n");
            for i = 1 : ndim
                Row             = self.proposalStartCorMat.Val(i,:);
                fprintf(outputUnit, formatVal + num2str(Row) + "\n");
            end
            Err.msg             = self.proposalStartCorMat.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit,  "\n" + "proposalStartStdVec" + "\n\n");
            for i = 1 : ndim
                fprintf(outputUnit, formatVal + num2str(self.proposalStartStdVec.Val(i)) + "\n");
            end
            Err.msg             = self.proposalStartStdVec.desc;
            if splashModeRequested, Err.note(); end

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