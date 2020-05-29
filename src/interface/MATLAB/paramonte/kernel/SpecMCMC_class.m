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

            Decoration          = Decoration_class([],[],[],[]);
            formatVal           = Decoration.TAB + Decoration.TAB;
            
            Err                 = Err_class();
            Err.prefix          = prefix;
            Err.outputUnit      = outputUnit;
            Err.resetEnabled    = false;

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "chainSize\n\n");
            fprintf(outputUnit,  formatVal + num2str(self.chainSize.val) + "\n");
            Err.msg             = self.chainSize.desc;
            if splashModeRequested, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomStartPointDomainLowerLimitVec" + "\n\n");
            for i = 1 : length(self.randomStartPointDomainLowerLimitVec.Val)
                fprintf(outputUnit,  formatVal + num2str(self.randomStartPointDomainLowerLimitVec.Val(i))+ "\n");
            end
            Err.msg             = self.randomStartPointDomainLowerLimitVec.desc;
            if splashModeRequested == true, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomStartPointDomainUpperLimitVec" + "\n\n");
            for i = 1 : length(self.randomStartPointDomainUpperLimitVec.Val)
                fprintf(outputUnit, formatVal + self.randomStartPointDomainUpperLimitVec.Val(i) + "\n");
            end
            Err.msg             = self.randomStartPointDomainUpperLimitVec.desc;
            if splashModeRequested == true, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "startPointVec" + "\n\n");
            for i = 1 : length(self.startPointVec.Val)
                fprintf(outputUnit, formatVal + self.startPointVec.Val(i) + "\n");
            end
            Err.msg             = self.startPointVec.desc;
            if splashModeRequested == true, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "randomStartPointRequested" + "\n\n");
            fprintf(outputUnit, formatVal + self.randomStartPointRequested.val + "\n");
            Err.msg             = self.randomStartPointRequested.desc;
            if splashModeRequested == true, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "sampleRefinementCount" + "\n\n");
            fprintf(outputUnit, formatVal + self.sampleRefinementCount.val + "\n");
            Err.msg             = self.sampleRefinementCount.desc;
            if splashModeRequested == true, Err.note(); end

            %-----------------------------------------------------------------------------------------------------------------------

            fprintf(outputUnit, "\n" + "sampleRefinementMethod" + "\n\n");
            fprintf(outputUnit, formatVal + self.sampleRefinementMethod.val + "\n");
            Err.msg             = self.sampleRefinementMethod.desc;
            if splashModeRequested == true, Err.note(); end
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