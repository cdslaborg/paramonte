classdef SpecBase_MaxNumDomainCheckToWarn_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_MaxNumDomainCheckToWarn_mod"
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

        function self = SpecBase_MaxNumDomainCheckToWarn_class()
            self.def    = 1000;
            self.desc   = "maxNumDomainCheckToWarn is an integer number beyond which the user will be warned about the newly-proposed points "          ...
                        + "being excessively proposed outside the domain of the objective function. For every maxNumDomainCheckToWarn "                 ...
                        + "consecutively-proposed new points that fall outside the domain of the objective function, the user will be warned until "    ...
                        + "maxNumDomainCheckToWarn = maxNumDomainCheckToStop, in which case the sampler returns a fatal error and the program stops "   ...
                        + "globally. The counter for this warning message is reset after a proposal sample from within the domain of the "              ...
                        + "objective function is obtained. The default value is " + num2str(self.def) + "."                                             ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, maxNumDomainCheckToWarn)  % set => setMaxNumDomainCheckToWarn
            self.val = maxNumDomainCheckToWarn;
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
                Err.msg         = Err.msg                                                                                                   ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                    ...
                                + "The input value for variable maxNumDomainCheckToWarn must be a positive integer. If you are not sure "   ...
                                + "about the appropriate value for this variable, simply drop it from the input. "                          ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline                 ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end