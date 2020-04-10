classdef SpecBase_OutputColumnWidth_class < handle

    properties (Constant)
        OUTPUT_COL_WIDTH        = 63
        CLASS_NAME              = "@SpecBase_OutputColumnWidth_class"
    end

    properties
        OutputColumnWidth_type  = []
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

        function self = SpecBase_OutputColumnWidth_class(methodName)
            self.def    = 0;
            self.desc   = "The variable outputColumnWidth is a non-negative integer number that determines the width of "                           ...
                        + "the data columns in " + methodName + " formatted output files that have tabular structure. "                             ...
                        + "If it is set to zero, " + methodName + " will ensure to set the width of each output element "                           ...
                        + "to the minimum possible width without losing the requested output precision. In other words, "                           ...
                        + "setting outputColumnWidth = 0 will result in the smallest-size for the formatted output files that are in ASCII format. "...
                        + "The default value is " + num2str(self.def) + "."                                                                         ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, outputColumnWidth)
            self.val = outputColumnWidth;
            if isempty(self.val)
                self.val = self.def;
            end
            self.str = num2str(self.val);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, outputRealPrecision)
            FUNCTION_NAME = "@checkForSanity()";
            if isempty(outputRealPrecision), outputRealPrecision = 0; end
            if self.val < 0
                Err.occurred    = true;
                Err.msg         = Errmsg                                                                                                    ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                    ...
                                + "The input value for variable outputColumnWidth must be a non-negative integer. "                         ...
                                + "If you are not sure about the appropriate value for this variable, simply drop it from the input. "      ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline                 ...
                                ;
            elseif (self.val > 0) && (self.val < (outputRealPrecision + 7))
                Err.occurred    = true;
                Err.msg         = Err.msg...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                    ...
                                + "The input value for variable outputColumnWidth must be outputColumnWidth > outputRealPrecision + 7. "    ...
                                + "If you are not sure about the appropriate value for this variable, either set it to zero on input, "     ...
                                + "or simply drop it from the input. " + methodName                                                         ...
                                + " will automatically assign an appropriate value to it." + newline + newline                              ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end