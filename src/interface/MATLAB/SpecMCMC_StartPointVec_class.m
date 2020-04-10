classdef SpecMCMC_StartPointVec_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecMCMC_StartPointVec_mod"
    end

    properties
        Val         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecMCMC_StartPointVec_class()
            self.desc   = "StartPoint is a 64bit real-valued vector of length ndim (the dimension of the domain of the input objective function)."  ...
                        + "For every element of StartPoint that is not provided as input, the default value will be the center of the domain of "   ...
                        + "StartPoint as specified by RandomStartPointDomain input variable. If the input variable RandomStartPointRequested=TRUE " ...
                        + "(or true or t, all case-insensitive), then the missing elements of StartPoint will be initialized to values drawn "      ...
                        + "randomly from within the corresponding ranges specified by the input variable RandomStartPointDomain."                   ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, startPointVec, randomStartPointDomainLowerLimitVec, randomStartPointDomainUpperLimitVec, randomStartPointRequested)
            if isempty(startPointVec)
                if randomStartPointRequested
                    self.Val = randomStartPointDomainLowerLimitVec + rand * (randomStartPointDomainUpperLimitVec - randomStartPointDomainLowerLimitVec);
                else
                    self.Val = 0.5 * (randomStartPointDomainLowerLimitVec + randomStartPointDomainUpperLimitVec);
                end
            else
                self.Val = startPointVec;
            end
            
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, randomStartPointDomainLowerLimitVec, randomStartPointDomainUpperLimitVec)
            FUNCTION_NAME = "@checkForSanity()";
            for i = 1 : length(self.Val)
                if (self.Val(i) < randomStartPointDomainLowerLimitVec(i)) || (self.Val(i) > randomStartPointDomainUpperLimitVec(i))
                    Err.occurred    = true;
                    Err.msg         = Err.msg...
                                    + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "...
                                    + "The input requested value for the component " + num2str(i) + " of the vector StartPoint ("   ...
                                    + num2str(self.Val(i)) + ") must be within the range of the sampling Domain defined "           ...
                                    + "in the program: ("                                                                           ...
                                    + num2str(randomStartPointDomainLowerLimitVec(i)) + ","                                         ...
                                    + num2str(randomStartPointDomainUpperLimitVec(i)) + "). If you don't "                          ...
                                    + "know an appropriate value for StartPoint, drop it from the input list. "                     ...
                                    + methodName + " will automatically assign an appropriate value to it." + newline + newline     ...
                                    ;
                end
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end