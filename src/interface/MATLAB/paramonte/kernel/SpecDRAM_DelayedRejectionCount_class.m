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

classdef SpecDRAM_DelayedRejectionCount_class < handle

    properties (Constant)
        CLASS_NAME                  = "SpecDRAM_DelayedRejectionCountclass"
        MAX_DELAYED_REJECTION_COUNT = 1000
        MIN_DELAYED_REJECTION_COUNT = 0
    end

    properties
        val                         = []
        def                         = []
        desc                        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_DelayedRejectionCount_class(methodName)
            self.def    = 0;
            self.desc   = num2str(self.MIN_DELAYED_REJECTION_COUNT) + " <= delayedRejectionCount <= " + num2str(self.MAX_DELAYED_REJECTION_COUNT)       ...
                        + " is an integer that represents the total number of stages for which rejections of new proposals "                            ...
                        + "will be tolerated by " + methodName + " before going back to the previously accepted point (state). "                        ...
                        + "Possible values are:" + newline + newline                                                                                    ...
                        + "    delayedRejectionCount = 0" + newline + newline                                                                           ...
                        + "            indicating no deployment of the delayed rejection algorithm." + newline + newline                                ...
                        + "    delayedRejectionCount > 0" + newline + newline                                                                           ...
                        + "            which implies a maximum delayedRejectionCount number of rejections will be tolerated." + newline + newline       ...
                        + "For example, delayedRejectionCount = 1, means that at any point during the sampling, if a proposal is rejected, "            ...
                        + methodName + " will not go back to the last sampled state. Instead, it will continue to propose a new from the current "      ...
                        + "rejected state. If the new state is again rejected based on the rules of " + methodName + ", then the algorithm will not "   ...
                        + "tolerate further rejections, because the maximum number of rejections to be tolerated has been set by the user to be "       ...
                        + "elayedRejectionCount = 1. The algorithm then goes back to the original last-accepted state and will begin proposing "        ...
                        + "new states from that location. "                                                                                             ...
                        + "The default value is delayedRejectionCount = " + num2str(self.def) + "."                                                     ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, delayedRejectionCount)
            if isempty(delayedRejectionCount)
                self.val = self.def;
            else
                self.val = delayedRejectionCount;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < self.MIN_DELAYED_REJECTION_COUNT
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                                ...
                                + "The input requested value for delayedRejectionCount (" + num2str(self.val) + ") "                                    ...
                                + "can not be negative. If you are not sure of the appropriate value for delayedRejectionCount, drop it "               ...
                                + "from the input list. " + methodName + " will automatically assign an appropriate value to it." + newline + newline   ...
                                ;
            elseif self.val > self.MAX_DELAYED_REJECTION_COUNT
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                                ...
                                + "The input requested value for delayedRejectionCount (" + num2str(self.val) + ") "                                    ...
                                + "can not be > " + num2str(self.MAX_DELAYED_REJECTION_COUNT)                                                           ...
                                + ". If you are not sure of the appropriate value for delayedRejectionCount, drop it "                                  ...
                                + "from the input list. " + methodName + " will automatically assign an appropriate value to it." + newline + newline   ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end