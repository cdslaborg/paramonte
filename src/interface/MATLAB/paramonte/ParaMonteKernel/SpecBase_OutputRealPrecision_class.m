classdef SpecBase_OutputRealPrecision_class < handle

    properties (Constant)
        CLASS_NAME              = "@SpecBase_OutputRealPrecision_class"
        OUTPUT_REAL_PRECISION   = 8
    end

    properties
        val                     = []
        def                     = []
        str                     = []
        desc                    = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_OutputRealPrecision_class(methodName)
            self.def    = self.OUTPUT_REAL_PRECISION;
            self.desc   = "The variable outputRealPrecision is a 32-bit integer number that determines the precision - that is, the number of "         ...
                        + "significant digits - of the real numbers in the output files of " + methodName + ". Any positive integer is acceptable "     ...
                        + "as the input value of outputRealPrecision. However, any digits of the output real numbers beyond the accuracy of 64-bit "    ...
                        + "real numbers (approximately 16 digits of significance) will be meaningless and random. Set this variable to 16 (or larger) " ...
                        + "if full reproducibility of the simulation is needed in the future. But keep in mind that larger precisions will result in "  ...
                        + "larger-size output files. This variable is ignored for binary output (if any occurs during the simulation). "                ...
                        + "The default value is " + num2str(self.def) + "."                                                                             ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, outputRealPrecision)
            self.val = outputRealPrecision;
            if isempty(self.val)
                self.val = self.def;
            end
            self.str = num2str(self.val);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME   = "@checkForSanity()";
            if self.val < 1
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                   ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                    ...
                                + "The input value for variable outputRealPrecision must be a positive integer < 16. If you are not sure "  ...
                                + "about the appropriate value for this variable, simply drop it from the input. "                          ...
                                + methodName + " will automatically assign an appropriate value to it." + newline +newline                  ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end