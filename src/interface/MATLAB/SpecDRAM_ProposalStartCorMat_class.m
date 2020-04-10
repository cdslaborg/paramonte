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