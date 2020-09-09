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

classdef SpecDRAM_ProposalStartCorMat_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_ProposalStartCorMat_class";
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

        function self = SpecDRAM_ProposalStartCorMat_class(nd, methodName)
            self.Def    = eye(nd,nd);
            self.desc   = "ProposalStartCorMat is a real-valued positive-definite matrix of size (ndim,ndim), where ndim is the dimension of the "              ...
                        + "sampling space. It serves as the best-guess starting correlation matrix of the multivariate Normal proposal distribution used by "   ...
                        + methodName + ". "                                                                                                                     ...
                        + "It is used (along with the input vector ProposalStartStdVec) to construct the covariance matrix of the proposal distribution "       ...
                        + "when the input covariance matrix is missing in the input list of variables. "                                                        ...
                        + "If the covariance matrix is given as input to " + methodName + ", any input values for "                                             ...
                        + "ProposalStartCorMat, as well as StartStdVec, will be automatically ignored by " + methodName + ". "                                  ...
                        + "As input to " + methodName + ", ProposalStartCorMat along with StdVec are especially useful when obtaining a best-guess "            ...
                        + "covariance matrix is not trivial. The default matrix value is simply the Identity Matrix of size ndim-by-ndim."                      ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, ProposalStartCorMat)
            if isempty(ProposalStartCorMat)
                self.Val = self.Def;
            else
                self.Val = ProposalStartCorMat;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, nd)
            FUNCTION_NAME   = "@checkForSanity()";
            
            [~,dummy]       = chol(self.Val(1:nd,1:nd));
            isPosDef        = ~dummy;
            if ~isPosDef
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred." + newline                       ...
                                + "The input requested ProposalStartCorMat for the Gaussian Proposal of " + methodName  ...
                                + " is not a positive-definite matrix." + newline + newline                             ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end