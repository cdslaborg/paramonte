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

classdef SpecBase_TargetAcceptanceRate_class < handle

    properties (Constant)
        CLASS_NAME          = "@SpecBase_TargetAcceptanceRate_class"
    end

    properties
        scalingRequested    = []
        val                 = []
        desc                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        function self = SpecBase_TargetAcceptanceRate_class(methodName)
            self.scalingRequested   = true;
            self.desc               = "targetAcceptanceRate sets an optimal target for the ratio of the number of accepted objective function calls to the "    ...
                                    + "total number of function calls by " + methodName + ". By default, it is a real number between 0 and 1. "                 ...
                                    + "If provided by the user, " + methodName + " will attempt (but not guarantee) to bring the average acceptance ratio "     ...
                                    + "of the sampler as close to the user-provided target ratio as possible. The success of " + methodName                     ...
                                    + " in keeping the average acceptance ratio close to the requested target value depends heavily on:" + newline              ...
                                    + "    1) the value of adaptiveUpdatePeriod; the larger, the easier." + newline                                             ...
                                    + "    2) the value of adaptiveUpdateCount; the larger, the easier." + newline                                              ...
                                    + "Note that the acceptance ratio adjustments will only occur every adaptiveUpdatePeriod sampling steps for a total "       ...
                                    + "number of adaptiveUpdateCount. "                                                                                         ...
                                    + "There is no default value for targetAcceptanceRate, as the acceptance ratio is not directly adjusted during sampling."   ...
                                    ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, targetAcceptanceRate)
            if isempty(targetAcceptanceRate)
                self.scalingRequested = false;
            else
                self.val = targetAcceptanceRate;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err)
            FUNCTION_NAME = "@checkForSanity()";
            if ~self.scalingRequested, return; end
            if self.val <= 0
                Err.occurred    = true;
                Err.msg         = Err.msg...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                    ...
                                + "The target acceptance ratio targetAcceptanceRate (" + num2str(self.val)  ...
                                + ") cannot be less than or equal to 0." + newline + newline                ...
                                ;
            end
            if self.val >= 1
                Err.occurred    = true;
                Err.msg         = Err.msg...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                    ...
                                + "The target acceptance ratio targetAcceptanceRate (" + num2str(self.val)  ...
                                + ") cannot be larger than or equal to 1." + newline + newline              ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end