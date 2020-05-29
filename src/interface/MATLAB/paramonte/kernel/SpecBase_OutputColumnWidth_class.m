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