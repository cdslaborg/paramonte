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
        overwriteRequested      = []
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
            self.overwriteRequested         = SpecBase_OverwriteRequested_class         (methodName);
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
                                    , overwriteRequested        ...
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

            self.overwriteRequested.        set(overwriteRequested);
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

            formatVal           = Decoration_class.TAB + Decoration_class.TAB;
            
            Err                 = Err_class();
            Err.prefix          = prefix;
            Err.outputUnit      = outputUnit;
            Err.resetEnabled    = false;

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "description" + "\n\n");
            Err.fullprefix      = formatVal;
            Err.msg             = self.description.val;
            Err.informUser();
            Err.msg             = self.description.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "silentModeRequested" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.silentModeRequested.val) + "\n");
            Err.msg             = self.silentModeRequested.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "domainLowerLimitVec" + "\n\n");
            for i = 1 : length(self.domainLowerLimitVec.Val)
                fprintf(outputUnit, formatVal + num2str(self.domainLowerLimitVec.Val(i)) + "\n");
            end
            Err.msg             = self.domainLowerLimitVec.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "domainUpperLimitVec" + "\n\n");
            for i = 1 : length(self.domainUpperLimitVec.Val)
                fprintf(outputUnit, formatVal + num2str(self.domainUpperLimitVec.Val(i)) + "\n");
            end
            Err.msg             = self.domainUpperLimitVec.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "variableNameList" + "\n\n");
            for i = 1: length(self.variableNameList.Val)
                fprintf(outputUnit, formatVal + char(self.variableNameList.Val{i}) + "\n");
            end
            Err.msg             = self.variableNameList.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "parallelizationModel" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.parallelizationModel.val) + "\n");
            Err.msg             = self.parallelizationModel.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputFileName" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(strrep(self.outputFileName.modified, '\', '\\')) + "\n");
            Err.msg             = self.outputFileName.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "overwriteRequested" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.overwriteRequested.val) + "\n");
            Err.msg             = self.overwriteRequested.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "targetAcceptanceRate" + "\n\n");
            if isempty(self.targetAcceptanceRate.val)
                fprintf(outputUnit, formatVal + "UNDEFINED" + "\n");
            else
                fprintf(outputUnit, formatVal + num2str(self.targetAcceptanceRate.val) + "\n");
            end
            Err.msg             = self.targetAcceptanceRate.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "sampleSize" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.sampleSize.val) + "\n");
            Err.msg             = self.sampleSize.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomSeed" + "\n\n");
            fprintf(outputUnit, formatVal + "User-Requested Value" + "\n");
            if isempty(self.randomSeed.seed)
                fprintf(outputUnit, formatVal + "Not provided by the user. The seed for each processor will be appropriately assigned." + "\n");
            else
                fprintf(outputUnit, formatVal + num2str(self.randomSeed.seed) + "\n");
            end
            Err.msg             = self.randomSeed.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputColumnWidth" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.outputColumnWidth.val) + "\n");
            Err.msg             = self.outputColumnWidth.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputDelimiter" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.outputDelimiter.val) + "\n");
            Err.msg             = self.outputDelimiter.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "outputRealPrecision" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.outputRealPrecision.val) + "\n");
            Err.msg             = self.outputRealPrecision.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "chainFileFormat" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.chainFileFormat.val) + "\n");
            Err.msg             = self.chainFileFormat.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "restartFileFormat" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.restartFileFormat.val) + "\n");
            Err.msg             = self.restartFileFormat.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "progressReportPeriod" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.progressReportPeriod.val) + "\n");
            Err.msg             = self.progressReportPeriod.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "maxNumDomainCheckToWarn" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.maxNumDomainCheckToWarn.val) + "\n");
            Err.msg             = self.maxNumDomainCheckToWarn.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "maxNumDomainCheckToStop" + "\n\n");
            fprintf(outputUnit, formatVal + num2str(self.maxNumDomainCheckToStop.val) + "\n");
            Err.msg             = self.maxNumDomainCheckToStop.desc;
            if self.silentModeRequested.isFalse, Err.note(); end

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