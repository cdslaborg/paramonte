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