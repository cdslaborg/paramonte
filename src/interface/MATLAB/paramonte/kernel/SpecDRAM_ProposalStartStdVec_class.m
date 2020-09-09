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

classdef SpecDRAM_ProposalStartStdVec_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_ProposalStartStdVec_class"
    end

    properties
        Val         = []
        Def         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_ProposalStartStdVec_class(nd, methodName)
            self.Def    = 0;
            for i = 1 : nd
                self.Def(i) = 1;
            end
            self.desc   = "ProposalStartStdVec is a real-valued positive vector of length ndim, where ndim is the dimension of the sampling space. "        ...
                        + "It serves as the best-guess starting Standard Deviation of each of the components of the proposal distribution. "                ...
                        + "If the initial covariance matrix (ProposalStartCovMat) is missing as an input variable to "                                      ...
                        + methodName + ", then ProposalStartStdVec (along with the input variable ProposalStartCorMat) will be used to construct "          ...
                        + "the initial covariance matrix of the proposal distribution of the MCMC sampler. "                                                ...
                        + "However, if ProposalStartCorMat is present as an input argument to " + methodName + ", then the input ProposalStartStdVec "      ...
                        + "(along with the input ProposalStartCovMat) will be completely ignored and the input value for ProposalStartCovMat will be used " ...
                        + "to initialize the initial covariance matrix of the proposal distribution for " + methodName + ". "                               ...
                        + "The default value of ProposalStartStdVec is a vector of unit values (i.e., ones) of length ndim."                                ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, ProposalStartStdVec)
            if isempty(ProposalStartStdVec)
                self.Val = self.Def;
            else
                self.Val = ProposalStartStdVec;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, nd)
            FUNCTION_NAME = "@checkForSanity()";
            for i = 1 : nd
                if self.Val(i) <= 0
                    Err.occurred    = true;
                    Err.msg         = Err.msg                                                                                   ...
                                    + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                    ...
                                    + "The input requested value (" + num2str(self.Val(i)) + ") for the component "             ...
                                    + num2str(i) + " of the variable ProposalStartStdVec for the Proposal distribution of "     ...
                                    + methodName + " must be a positive real number." + newline + newline                       ...
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