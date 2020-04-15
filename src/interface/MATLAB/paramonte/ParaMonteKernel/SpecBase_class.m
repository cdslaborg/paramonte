classdef SpecBase_class < handle

    properties (Constant)
        CLASS_NAME = "@SpecBase_class"
    end

    properties
        sampleSize              = []
        randomSeed              = []
        description             = []
        outputFileName          = []
        outputDelimiter         = []
        chainFileFormat         = []
        variableNameList        = []
        domainLowerLimitVec     = []
        domainUpperLimitVec     = []
        restartFileFormat       = []
        outputColumnWidth       = []
        outputRealPrecision     = []
        silentModeRequested     = []
        parallelizationModel    = []
        progressReportPeriod    = []
        targetAcceptanceRate    = []
        maxNumDomainCheckToWarn = []
        maxNumDomainCheckToStop = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_class(nd,methodName)
            self.sampleSize                 = SpecBase_SampleSize_class                 (methodName);
            self.randomSeed                 = SpecBase_RandomSeed_class                 (methodName);
            self.description                = SpecBase_Description_class                (methodName);
            self.outputFileName             = SpecBase_OutputFileName_class             (methodName);
            self.outputDelimiter            = SpecBase_OutputDelimiter_class            (methodName);
            self.chainFileFormat            = SpecBase_ChainFileFormat_class            (methodName);
            self.variableNameList           = SpecBase_VariableNameList_class           (methodName, nd);
            self.restartFileFormat          = SpecBase_RestartFileFormat_class          (methodName);
            self.outputColumnWidth          = SpecBase_OutputColumnWidth_class          (methodName);
            self.outputRealPrecision        = SpecBase_OutputRealPrecision_class        (methodName);
            self.silentModeRequested        = SpecBase_SilentModeRequested_class        (methodName);
            self.domainLowerLimitVec        = SpecBase_DomainLowerLimitVec_class        (methodName);
            self.domainUpperLimitVec        = SpecBase_DomainUpperLimitVec_class        (methodName);
            self.parallelizationModel       = SpecBase_ParallelizationModel_class       (methodName);
            self.progressReportPeriod       = SpecBase_ProgressReportPeriod_class       ( );
            self.targetAcceptanceRate       = SpecBase_TargetAcceptanceRate_class       (methodName);
            self.maxNumDomainCheckToWarn    = SpecBase_MaxNumDomainCheckToWarn_class    ( );
            self.maxNumDomainCheckToStop    = SpecBase_MaxNumDomainCheckToStop_class    ( );
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function setFromInputArgs   ( self                      ...
                                    , ndim                      ...
                                    , Err                       ...
                                    , sampleSize                ...
                                    , randomSeed                ...
                                    , description               ...
                                    , outputFileName            ...
                                    , outputDelimiter           ...
                                    , chainFileFormat           ...
                                    , variableNameList          ...
                                    , restartFileFormat         ...
                                    , outputColumnWidth         ...
                                    , outputRealPrecision       ...
                                    , silentModeRequested       ...
                                    , domainLowerLimitVec       ...
                                    , domainUpperLimitVec       ...
                                    , parallelizationModel      ...
                                    , progressReportPeriod      ...
                                    , targetAcceptanceRate      ...
                                    , maxNumDomainCheckToWarn   ...
                                    , maxNumDomainCheckToStop   ...
                                    )

            FUNCTION_NAME = self.CLASS_NAME + "@setFromInputArgs()";

            self.sampleSize.                set(sampleSize);
            self.randomSeed.                set(randomSeed);
            self.description.               set(description);
            self.outputFileName.            set(outputFileName);
            self.chainFileFormat.           set(chainFileFormat);
            self.variableNameList.          set(variableNameList);
            self.restartFileFormat.         set(restartFileFormat);
            self.domainLowerLimitVec.       set(ndim, domainLowerLimitVec);
            self.domainUpperLimitVec.       set(ndim, domainUpperLimitVec);

            % do not change the order with outputDelimiter
            self.outputColumnWidth.         set(outputColumnWidth);
            self.outputDelimiter.           set(outputDelimiter, self.outputColumnWidth.val);

            self.parallelizationModel.      set(parallelizationModel);
            self.outputRealPrecision.       set(outputRealPrecision);
            self.silentModeRequested.       set(silentModeRequested);
            self.progressReportPeriod.      set(progressReportPeriod);
            self.targetAcceptanceRate.      set(targetAcceptanceRate);
            self.maxNumDomainCheckToWarn.   set(maxNumDomainCheckToWarn);
            self.maxNumDomainCheckToStop.   set(maxNumDomainCheckToStop);

            if Err.occurred, Err.msg = FUNCTION_NAME + Err.msg; end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function reportValues(self, prefix, outputUnit)

            formatVal   = Decoration_class.TAB + Decoration_class.TAB;

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "description" + "\n\n");
            self.Err.msg = self.description.val;
            self.Err.prefix = Decoration_class.TAB + Decoration_class.TAB;
            self.Err.outputUnit = outputUnit;
            self.Err.wrapWidth = 125;
            self.Err.marginTop = 0;
            self.Err.marginBot = 0;
            self.Err.informUser();
            if self.silentModeRequested.isFalse, Err_class.note(self.description.desc, prefix, '', outputUnit, '', []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "silentModeRequested" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.silentModeRequested.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.silentModeRequested.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "domainLowerLimitVec" + "\n\n");
            for i = 1 : length(self.domainLowerLimitVec.Val)
                fprintf(outputUnit, formatVal + num2str(self.domainLowerLimitVec.Val(i)) + "\n");
            end
            if self.silentModeRequested.isFalse, Err_class.note(self.domainLowerLimitVec.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "domainUpperLimitVec" + "\n\n");
            for i = 1 : length(self.domainUpperLimitVec.Val)
                fprintf(outputUnit, formatVal + num2str(self.domainUpperLimitVec.Val(i)) + "\n");
            end
            if self.silentModeRequested.isFalse, Err_class.note(self.domainUpperLimitVec.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "variableNameList" + "\n\n");
            for i = 1: length(self.variableNameList.Val)
                fprintf(outputUnit, formatVal + char(self.variableNameList.Val{i}) + "\n");
            end
            if self.silentModeRequested.isFalse, Err_class.note(self.variableNameList.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "parallelizationModel" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.parallelizationModel.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.parallelizationModel.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputFileName" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(strrep(self.outputFileName.modified, '\', '\\')) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.outputFileName.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "targetAcceptanceRate" + "\n\n");
            if isempty(self.targetAcceptanceRate.val)
                fprintf(outputUnit, formatVal + "UNDEFINED" + "\n");
            else
                fprintf(outputUnit, formatVal + num2str(self.targetAcceptanceRate.val) + "\n");
            end
            if self.silentModeRequested.isFalse, Err_class.note(self.targetAcceptanceRate.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "sampleSize" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.sampleSize.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.sampleSize.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomSeed" + "\n\n");
            fprintf(outputUnit, formatVal + "User-Requested Value" + "\n");
            if isempty(self.randomSeed.seed)
                fprintf(outputUnit, formatVal + "Not provided by the user. The seed for each processor will be appropriately assigned." + "\n");
            else
                fprintf(outputUnit, formatVal + num2str(self.randomSeed.seed) + "\n");
            end
            if self.silentModeRequested.isFalse, Err_class.note(self.randomSeed.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputColumnWidth" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.outputColumnWidth.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.outputColumnWidth.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputDelimiter" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.outputDelimiter.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.outputDelimiter.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputRealPrecision" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.outputRealPrecision.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.outputRealPrecision.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "chainFileFormat" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.chainFileFormat.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.chainFileFormat.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "restartFileFormat" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.restartFileFormat.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.restartFileFormat.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "progressReportPeriod" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.progressReportPeriod.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.progressReportPeriod.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "maxNumDomainCheckToWarn" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.maxNumDomainCheckToWarn.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.maxNumDomainCheckToWarn.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "maxNumDomainCheckToStop" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.maxNumDomainCheckToStop.val) + "\n");
            if self.silentModeRequested.isFalse, Err_class.note(self.maxNumDomainCheckToStop.desc, prefix, '', outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, methodName, Err)
        
            Err = self.chainFileFormat.             checkForSanity(Err, methodName);
            Err = self.outputDelimiter.             checkForSanity(Err, methodName);
            Err = self.restartFileFormat.           checkForSanity(Err, methodName);
            Err = self.outputColumnWidth.           checkForSanity(Err, methodName, self.outputRealPrecision.val);
            Err = self.outputRealPrecision.         checkForSanity(Err, methodName);
            Err = self.domainLowerLimitVec.         checkForSanity(Err);
            Err = self.domainUpperLimitVec.         checkForSanity(Err,self.domainLowerLimitVec.Val);
            Err = self.parallelizationModel.        checkForSanity(Err, methodName);
            Err = self.progressReportPeriod.        checkForSanity(Err, methodName);
            Err = self.targetAcceptanceRate.        checkForSanity(Err);
            Err = self.maxNumDomainCheckToWarn.     checkForSanity(Err, methodName);
            Err = self.maxNumDomainCheckToStop.     checkForSanity(Err, methodName);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end