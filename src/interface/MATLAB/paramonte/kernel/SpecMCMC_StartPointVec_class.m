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

classdef SpecMCMC_StartPointVec_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecMCMC_StartPointVec_mod"
    end

    properties
        Val         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecMCMC_StartPointVec_class()
            self.desc   = "StartPoint is a 64bit real-valued vector of length ndim (the dimension of the domain of the input objective function)."  ...
                        + "For every element of StartPoint that is not provided as input, the default value will be the center of the domain of "   ...
                        + "StartPoint as specified by RandomStartPointDomain input variable. If the input variable RandomStartPointRequested=TRUE " ...
                        + "(or true or t, all case-insensitive), then the missing elements of StartPoint will be initialized to values drawn "      ...
                        + "randomly from within the corresponding ranges specified by the input variable RandomStartPointDomain."                   ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, startPointVec, randomStartPointDomainLowerLimitVec, randomStartPointDomainUpperLimitVec, randomStartPointRequested)
            if isempty(startPointVec)
                if randomStartPointRequested
                    self.Val = randomStartPointDomainLowerLimitVec + rand * (randomStartPointDomainUpperLimitVec - randomStartPointDomainLowerLimitVec);
                else
                    self.Val = 0.5 * (randomStartPointDomainLowerLimitVec + randomStartPointDomainUpperLimitVec);
                end
            else
                self.Val = startPointVec;
            end
            
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, randomStartPointDomainLowerLimitVec, randomStartPointDomainUpperLimitVec)
            FUNCTION_NAME = "@checkForSanity()";
            for i = 1 : length(self.Val)
                if (self.Val(i) < randomStartPointDomainLowerLimitVec(i)) || (self.Val(i) > randomStartPointDomainUpperLimitVec(i))
                    Err.occurred    = true;
                    Err.msg         = Err.msg...
                                    + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "...
                                    + "The input requested value for the component " + num2str(i) + " of the vector StartPoint ("   ...
                                    + num2str(self.Val(i)) + ") must be within the range of the sampling Domain defined "           ...
                                    + "in the program: ("                                                                           ...
                                    + num2str(randomStartPointDomainLowerLimitVec(i)) + ","                                         ...
                                    + num2str(randomStartPointDomainUpperLimitVec(i)) + "). If you don't "                          ...
                                    + "know an appropriate value for StartPoint, drop it from the input list. "                     ...
                                    + methodName + " will automatically assign an appropriate value to it." + newline + newline     ...
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