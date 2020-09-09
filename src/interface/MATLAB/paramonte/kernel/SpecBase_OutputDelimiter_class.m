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

classdef SpecBase_OutputDelimiter_class < handle

    properties (Constant)
        CLASS_NAME          = "@SpecBase_OutputDelimiter_class"
        MAX_DELIMITER_LEN   = 63
    end

    properties
        val                 = []
        def                 = []
        desc                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_OutputDelimiter_class(methodName)
            self.def    = ",";
            self.desc   = "outputDelimiter is a string variable, containing a sequence of one or more characters, that is used to specify the " ...
                        + "boundary between separate, independent information elements in the tabular output files of " + methodName + ". "     ...
                        + "The string value must be enclosed by either single or double quotation marks when provided as input. "               ...
                        + "To output in Comma-Separated-Values (CSV) format, set outputDelimiter = ','. If the input value is not provided, "   ...
                        + "the default delimiter '" + self.def + "' will be used when input outputColumnWidth = 0, and a single "               ...
                        + "space character, '" + self.def + "' will be used when input outputColumnWidth > 0. "                                 ...
                        + "The default value is '" + self.def + "'."                                                                            ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, outputDelimiter, outputColumnWidth)   % set => setOutputDelimiter
            if isempty(outputDelimiter)
                if outputColumnWidth == 0
                    self.val = self.def;
                else
                    self.val = " ";
                end
            else
                self.val = strtrim(outputDelimiter);
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME   = "@checkForSanity()";
            delimiter       = char(strtrim(self.val));
            delimiterLen    = length(delimiter);
            for i = 1 : delimiterLen
                if (isstrprop(delimiter(i),'digit')) || (delimiter(i) == ".") || (delimiter(i:i) == "-") || (delimiter(i:i) == "+")
                    Err.occurred = true;
                    break
                end
            end
            if Err.occurred
                Err.msg = Err.msg                                                                                                   ...
                        + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                    ...
                        + "The input value for variable outputDelimiter cannot contain any digits or the period symbol '.' or '-' " ...
                        + "or '+'. If you are unsure about the appropriate value for this variable, simply drop it from the input." ...
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