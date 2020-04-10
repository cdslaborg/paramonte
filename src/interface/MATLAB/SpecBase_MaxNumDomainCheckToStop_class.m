classdef SpecBase_MaxNumDomainCheckToStop_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_MaxNumDomainCheckToStop_class"
    end

    properties
        val         = []
        def         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_MaxNumDomainCheckToStop_class()
            self.def    = 10000;
            self.desc   = "maxNumDomainCheckToStop is an integer number beyond which the program will stop globally with a fatal error message "        ...
                        + "declaring that the maximum number of proposal-out-of-domain-bounds has reached. The counter for this global-stop request "   ...
                        + "is reset after a proposal point is accepted as a sample from within the domain of the objective function. "                  ...
                        + "The default value is " + num2str(self.def) + "."                                                                             ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, maxNumDomainCheckToStop)
            self.val = maxNumDomainCheckToStop;
            if isempty(self.val)
                self.val = self.def;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < 1
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                ...
                                + "The input value for variable maxNumDomainCheckToStop must be a positive integer. "                   ...
                                + "If you are not sure about the appropriate value for this variable, simply drop it from the input. "  ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline             ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
