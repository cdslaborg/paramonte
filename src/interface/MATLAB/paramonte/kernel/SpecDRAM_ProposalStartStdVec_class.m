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