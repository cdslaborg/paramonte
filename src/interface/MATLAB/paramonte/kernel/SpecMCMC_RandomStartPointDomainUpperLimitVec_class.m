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

classdef SpecMCMC_RandomStartPointDomainUpperLimitVec_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecMCMC_RandomStartPointDomainUpperLimitVec_class"
    end

    properties
        desc        = []
        Val         = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecMCMC_RandomStartPointDomainUpperLimitVec_class(methodName)
            self.desc   = "randomStartPointDomainLowerLimitVec represents the upper boundaries of the cubical domain from which the starting point(s) of "      ...
                        + "the MCMC chain(s) will be initialized randomly (only if requested via the input variable randomStartPointRequested. "                ...
                        + "This happens only when some or all of the elements of the input variable StartPoint are missing. "                                   ...
                        + "In such cases, every missing value of input StartPoint will be set to the center point between randomStartPointDomainLowerLimitVec " ...
                        + "and RandomStartPointDomainLowerLimitVec in the corresponding dimension. "                                                            ...
                        + "If RandomStartPointRequested=TRUE (or True, true, t, all case-insensitive), then the missing "                                       ...
                        + "elements of StartPoint will be initialized to values drawn randomly from within the corresponding ranges specified by "              ...
                        + "the input variable randomStartPointDomainLowerLimitVec. "                                                                            ...
                        + "As an input variable, randomStartPointDomainLowerLimitVec is an ndim-dimensional vector of 64-bit real numbers, "                    ...
                        + "where ndim is the number of variables of the objective function. It is also possible to assign only select values of "               ...
                        + "randomStartPointDomainLowerLimitVec and leave the rest of the components to be assigned the default value. "                         ...
                        + "This is POSSIBLE ONLY when randomStartPointDomainLowerLimitVec is defined inside the input file to " + methodName + ". "             ...
                        + "For example, having the following inside the input file, " + newline + newline                                                       ...
                        + "    randomStartPointDomainLowerLimitVec(3:5) = -100" + newline + newline                                                             ...
                        + "        will only set the upper limits of the third, fourth, and the fifth dimensions to -100, or," + newline + newline              ...
                        + "    randomStartPointDomainLowerLimitVec(1) = -100, randomStartPointDomainLowerLimitVec(2) = -1.e6 " + newline + newline              ...
                        + "        will set the upper limit on the first dimension to -100, and 1.e6 on the second dimension, or," + newline + newline          ...
                        + "    randomStartPointDomainLowerLimitVec = 3*-2.5e100" + newline + newline                                                            ...
                        + "        will only set the upper limits on the first, second, and the third dimensions to -2.5*10^100, while the rest of "            ...
                        +         "the upper limits for the missing dimensions will be automatically set to the default value." + newline + newline             ...
                        + "The default values for all elements of randomStartPointDomainLowerLimitVec are taken from the corresponding values in the input "    ...
                        + "variable DomainUpperLimitVec."                                                                                                       ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, randomStartPointDomainUpperLimitVec, DomainUpperLimitVec)
            if isempty(randomStartPointDomainUpperLimitVec)
                self.Val = DomainUpperLimitVec;
            else
                self.Val = randomStartPointDomainUpperLimitVec;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, RandomStartPointDomainLowerLimitVec, DomainUpperLimitVec)
            FUNCTION_NAME = "@checkForSanity()";
            for i = 1 : length(self.Val)
                if self.Val(i) > DomainUpperLimitVec(i)
                    Err.occurred    = true;
                    Err.msg         = Err.msg                                                                                               ...
                                    + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                ...
                                    + "The component " + num2str(i) + " of the variable randomStartPointDomainUpperLimitVec ("              ...
                                    + num2str(self.Val(i))                                                                                  ...
                                    + ") cannot be larger than the corresponding component of the variable "                                ...
                                    + "domainUpperLimitVec (" + num2str(DomainUpperLimitVec(i)) + "). If you don't know "                   ...
                                    + "an appropriate value to set for randomStartPointDomainUpperLimitVec, drop it from the input list. "  ...
                                    + methodName + " will automatically assign an appropriate value to it." + newline + newline             ...
                                    ;
                end
                if self.Val(i) <= RandomStartPointDomainLowerLimitVec(i)
                    Err.occurred    = true;
                    Err.msg         = Err.msg                                                                                                                       ...
                                    + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. The input upper limit value in the component " + num2str(i)              ...
                                    + " of the variable randomStartPointDomainUpperLimitVec cannot be smaller than or equal to the corresponding input "            ...
                                    + "lower limit value in randomStartPointDomainLowerLimitVec:\n"                                                                 ...
                                    + "    randomStartPointDomainLowerLimitVec(" + num2str(i) + ") = " + num2str(RandomStartPointDomainLowerLimitVec(i)) + newline  ...
                                    + "    randomStartPointDomainUpperLimitVec(" + num2str(i) + ") = " + num2str(self.Val(i)) + newline + newline                   ...
                                    ;
                end
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end