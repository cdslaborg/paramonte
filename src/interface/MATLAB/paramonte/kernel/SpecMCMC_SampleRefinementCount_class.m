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

classdef SpecMCMC_SampleRefinementCount_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecMCMC_SampleRefinementCount_class"
    end

    properties
        val         = []
        def         = []
        str         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        function self = SpecMCMC_SampleRefinementCount_class(methodName)
            self.def    = Constants.POSINF_IK;
            self.desc   = "When sampleSize < 0, the variable sampleRefinementCount is an integer that dictates the maximum number of times "                ...
                        + "the MCMC chain will be refined to remove the autocorrelation within the output MCMC sample. For example," + newline + newline    ...
                        + "    if sampleRefinementCount = 0,"+ newline + newline                                                                            ...
                        + "            no refinement of the output MCMC chain will be performed, the resulting MCMC sample will simply correspond "         ...
                        +             "to the full MCMC chain in verbose format (i.e., each sampled state has a weight of one)." + newline + newline        ...
                        + "    if sampleRefinementCount = 1," + newline + newline                                                                           ...
                        + "            the refinement of the output MCMC chain will be done only once if needed, and no more, "                             ...
                        +             "even though there may still exist some residual autocorrelation in the output MCMC sample. "                         ...
                        +             "In practice, only one refinement of the final output MCMC Chain should be enough to remove "                         ...
                        +             "the existing autocorrelations in the final output sample. Exceptions occur when the Integrated "                     ...
                        +             "Autocorrelation (IAC) of the output MCMC chain is comparable to or larger than the length of the chain. "            ...
                        +             "In such cases, neither the BatchMeans method nor any other method of IAC computation will be able to "               ...
                        +             "accurately compute the IAC. Consequently, the samples generated based on the computed IAC values will "              ...
                        +             "likely not be i.i.d. and will still be significantly autocorrelated. In such scenarios, more than "                  ...
                        +             "one refinement of the MCMC chain will be necessary. Very small sample size resulting from multiple "                 ...
                        +             "refinements of the sample could be a strong indication of the bad mixing of the MCMC chain and "                     ...
                        +             "the output chain may not contain true i.i.d. samples from the target objective function." + newline + newline        ...
                        + "    if sampleRefinementCount > 1," + newline + newline                                                                           ...
                        + "            the refinement of the output MCMC chain will be done for a maximum sampleRefinementCount number of times, "          ...
                        +             "even though there may still exist some residual autocorrelation in the final output MCMC sample." + newline + newline...
                        + "    if sampleRefinementCount >> 1 (e.g., comparable to or larger than the length of the MCMC chain)," + newline + newline        ...
                        + "            the refinement of the output MCMC chain will continue until the integrated autocorrelation of the resulting "        ...
                        +             "final sample is less than 2, virtually implying that an independent identically-distributed (i.i.d.) sample "        ...
                        +             "has finally been obtained." + newline + newline                                                                      ...
                        + "Note that to obtain i.i.d. samples from a multidimensional chain, " + methodName + " will use the maximum of "                   ...
                        + "Integrated Autocorrelation (IAC) among all dimensions of the chain to refine the chain. "                                        ...
                        + "Note that the value specified for sampleRefinementCount is used only when the variable sampleSize < 0, otherwise, "              ...
                        + "it will be ignored. The default value is sampleRefinementCount = " +  num2str(self.def) + "."                                    ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, sampleRefinementCount)
            self.val = sampleRefinementCount;
            if isempty(self.val)
                self.val = self.def;
            end
            self.str = num2str(self.val);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < 0 || floor(self.val) ~= self.val
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                ...
                                + "The input value for variable sampleRefinementCount must be a non-negative integer. "                 ...
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