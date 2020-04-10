classdef SpecMCMC_class < handle

    properties (Constant)
    end

    properties
        chainSize                               = []
        startPointVec                           = []
        sampleRefinementCount                   = []
        sampleRefinementMethod                  = []
        randomStartPointRequested               = []
        randomStartPointDomainLowerLimitVec     = []
        randomStartPointDomainUpperLimitVec     = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %********************************************************************************************************************************

        function self = SpecMCMC_class(methodName)

            self.chainSize                              = SpecMCMC_ChainSize_class(methodName);
            self.startPointVec                          = SpecMCMC_StartPointVec_class();
            self.sampleRefinementCount                  = SpecMCMC_SampleRefinementCount_class(methodName);
            self.sampleRefinementMethod                 = SpecMCMC_SampleRefinementMethod_class(methodName);
            self.randomStartPointRequested              = SpecMCMC_RandomStartPointRequested_class(methodName);
            self.randomStartPointDomainLowerLimitVec    = SpecMCMC_RandomStartPointDomainLowerLimitVec_class(methodName);
            self.randomStartPointDomainUpperLimitVec    = SpecMCMC_RandomStartPointDomainUpperLimitVec_class(methodName);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function setFromInputArgs ( self, domainLowerLimitVec, domainUpperLimitVec  ...
                                    ...ParaMCMC variables
                                    , chainsize                                     ...
                                    , startPointVec                                 ...
                                    , sampleRefinementCount                         ...
                                    , sampleRefinementMethod                        ...
                                    , randomStartPointRequested                     ...
                                    , randomStartPointDomainLowerLimitVec           ...
                                    , randomStartPointDomainUpperLimitVec           ...
                                    )
            self.chainSize.                             set ( chainsize);
            self.sampleRefinementMethod.                set ( sampleRefinementMethod);
            self.randomStartPointRequested.             set ( randomStartPointRequested);
            self.sampleRefinementCount.                 set ( sampleRefinementCount);
            self.randomStartPointDomainLowerLimitVec.   set ( randomStartPointDomainLowerLimitVec, domainLowerLimitVec);
            self.randomStartPointDomainUpperLimitVec.   set ( randomStartPointDomainUpperLimitVec, domainUpperLimitVec);
            self.startPointVec.                            set ( startPointVec, self.randomStartPointDomainLowerLimitVec.Val    ...
                                                            , self.randomStartPointDomainUpperLimitVec.Val                      ...
                                                            , self.randomStartPointRequested.val                                ...
                                                            );
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function reportValues(self, prefix, outputUnit, splashModeRequested)

            Decoration = Decoration_class([],[],[],[]);
            formatVal = Decoration.TAB + Decoration.TAB;

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "chainSize\n\n");
            fprintf(outputUnit,  formatVal + num2str(self.chainSize.val) + "\n");
            if splashModeRequested, Err_class.note(self.chainSize.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomStartPointDomainLowerLimitVec" + "\n\n");
            for i = 1 : length(self.randomStartPointDomainLowerLimitVec.Val)
                fprintf(outputUnit,  formatVal + num2str(self.randomStartPointDomainLowerLimitVec.Val(i))+ "\n");
            end
            if splashModeRequested == true, Err_class.note(self.randomStartPointDomainLowerLimitVec.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomStartPointDomainUpperLimitVec" + "\n\n");
            for i = 1 : length(self.randomStartPointDomainUpperLimitVec.Val)
                fprintf(outputUnit, formatVal + self.randomStartPointDomainUpperLimitVec.Val(i) + "\n");
            end
            if splashModeRequested == true, Err_class.note(self.randomStartPointDomainUpperLimitVec.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "startPointVec" + "\n\n");
            for i = 1 : length(self.startPointVec.Val)
                fprintf(outputUnit, formatVal + self.startPointVec.Val(i) + "\n");
            end
            if splashModeRequested == true, Err_class.note(self.startPointVec.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomStartPointRequested" + "\n\n");
            fprintf(outputUnit, formatVal + self.randomStartPointRequested.val + "\n");
            if splashModeRequested == true, Err_class.note(self.randomStartPointRequested.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "sampleRefinementCount" + "\n\n");
            fprintf(outputUnit, formatVal + self.sampleRefinementCount.val + "\n");
            if splashModeRequested == true, Err_class.note(self.sampleRefinementCount.desc, prefix, "\n", outputUnit, [], []); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "sampleRefinementMethod" + "\n\n");
            fprintf(outputUnit, formatVal + self.sampleRefinementMethod.val + "\n");
            if splashModeRequested == true, Err_class.note(self.sampleRefinementMethod.desc, prefix, "\n", outputUnit, [], []); end
            %-----------------------------------------------------------------------------------------------------------------------

        end % function reportValues

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, nd, domainLowerLimitVec, domainUpperLimitVec)
            self.chainSize.                             checkForSanity  ( Err, methodName, nd);
            self.sampleRefinementCount.                 checkForSanity  ( Err, methodName);
            self.sampleRefinementMethod.                checkForSanity  ( Err, methodName);
            self.randomStartPointDomainLowerLimitVec.   checkForSanity  ( Err, methodName, domainLowerLimitVec);
            self.randomStartPointDomainUpperLimitVec.   checkForSanity  ( Err, methodName, self.randomStartPointDomainLowerLimitVec.Val, domainUpperLimitVec);
            self.startPointVec.                         checkForSanity  ( Err, methodName                               ...
                                                                        , self.randomStartPointDomainLowerLimitVec.Val  ...
                                                                        , self.randomStartPointDomainUpperLimitVec.Val  ...
                                                                        );
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end