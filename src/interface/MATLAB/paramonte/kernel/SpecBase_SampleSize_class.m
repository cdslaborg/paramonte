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

classdef SpecBase_SampleSize_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_SampleSize_class"
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

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_SampleSize_class(methodName)
            self.def    = -1;
            self.desc   = "The variable sampleSize is an integer that dictates the number of (hopefully, independent and identically distributed "              ...
                        + "[i.i.d.]) samples to be drawn from the user-provided objective function. Three ranges of values are possible:" + newline + newline   ...
                        + "    sampleSize < 0:" + newline + newline                                                                                             ...
                        + "            Then, the absolute value of sampleSize dictates the sample size in units of the effective sample size. "                 ...
                        + "The effective sample is by definition i.i.d., and free from duplicates. The effective sample size "                                  ...
                        + "is determined by " + methodName + " automatically toward the end of the simulation."                                                 ...
                        + "            For example:" + newline + newline                                                                                        ...
                        + "                    sampleSize = -1 yields the effective i.i.d. sample drawn from the objective function." + newline + newline       ...
                        + "                    sampleSize = -2 yields a (potentially non-i.i.d.) sample twice as big as the effective "                         ...
                        + "sample." + newline + newline                                                                                                         ...
                        + "    sampleSize > 0:" + newline + newline                                                                                             ...
                        + "            Then, the sample size is assumed to be in units of the number of points to be sampled. "                                 ...
                        + "If sampleSize turns out to be less than effectiveSampleSize, the resulting sample will be i.i.d.. "                                  ...
                        + "If sampleSize turns out to be larger than effectiveSampleSize, the resulting sample will be "                                        ...
                        + "potentially non-i.i.d.. The larger the difference, the more non-i.i.d. the resulting sample will be." + newline                      ...
                        + "            For example:" + newline + newline                                                                                        ...
                        + "                    sampleSize = 1000 yields a 1000-points sample from the objective function." + newline + newline                  ...
                        + "    sampleSize = 0:" + newline + newline                                                                                             ...
                        + "            in which case, no sample file will be generated." + newline + newline                                                    ...
                        + "Default value is sampleSize = " + num2str(self.def) + "."                                                                            ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, sampleSize)
            if isempty(sampleSize)
                self.val = self.def;
            else
                self.val = sampleSize;
            end
            self.str = num2str(self.val);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self,Err,methodName)
            if self.val < 1
                Err.occurred    = true;
                Err.msg         = Err.msg + self.CLASS_NAME + "@checkForSanity(): Error occurred. "                                     ...
                                + "The input value for variable sampleSize must be a positive integer. If you are not sure about the "  ...
                                + "appropriate value for this variable, simply drop it from the input. " + methodName                   ...
                                + " will automatically assign an appropriate value to it." + newline + newline                          ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end