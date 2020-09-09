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

classdef SpecDRAM_BurninAdaptationMeasure_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_BurninAdaptationMeasure_class"
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

        function self = SpecDRAM_BurninAdaptationMeasure_class(methodName)
            self.def    = 1;
            self.desc   = "burninAdaptationMeasure is a 64-bit real number between 0 and 1, representing the adaptation measure threshold below "               ...
                        + "which the simulated Markov chain will be used to generate the output " + methodName + " sample. In other words, any "                ...
                        + "point in the output Markov Chain that has been sampled during significant adaptation of the proposal distribution (as "              ...
                        + "determined by burninAdaptationMeasure) will not be included in the construction of the final " + methodName + " output "             ...
                        + "sample. This is to ensure that the generation of the output sample will be based on the part of the simulated chain that "           ...
                        + "is practically guaranteed to be Markovian and ergodic. If this variable is set to 0, then the output sample will be "                ...
                        + "generated from the part of the chain where no proposal adaptation has occurred. This non-adaptive or "                               ...
                        + "minimally-adaptive part of the chain may not even exist if the total adaptation period of the simulation "                           ...
                        + "(as determined by adaptiveUpdateCount and adaptiveUpdatePeriod input variables) is longer than the total length "                    ...
                        + "of the output MCMC chain. In such cases, the resulting output sample may have a zero size. In general, when "                        ...
                        + "good mixing occurs (e.g., when the input variable chainSize is very large) any specific value of burninAdaptationMeasure "           ...
                        + "becomes practically irrelevant. "                                                                                                    ...
                        + "The default value for burninAdaptationMeasure is " + num2str(self.def) + ", implying that the entire "                               ...
                        + "chain (with the exclusion of an initial automatically-determined burnin period) will be used to generate the final output sample."   ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, burninAdaptationMeasure)
            if isempty(burninAdaptationMeasure)
                self.val = self.def;
            else
                self.val = burninAdaptationMeasure;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < 0
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                                ...
                                + "The input variable burninAdaptationMeasure (" + num2str(self.val)                                                    ...
                                + ") cannot be less than 0. If you are not sure of the appropriate value for burninAdaptationMeasure, drop it "         ...
                                + "from the input list. " + methodName + " will automatically assign an appropriate value to it" + newline + newline    ...
                                ;
            end
            if self.val > 1
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                                ...
                                + "The input variable burninAdaptationMeasure (" + num2str(self.val)                                                    ...
                                + ") cannot be larger than 1. If you are not sure of the appropriate value for burninAdaptationMeasure, drop it "       ...
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