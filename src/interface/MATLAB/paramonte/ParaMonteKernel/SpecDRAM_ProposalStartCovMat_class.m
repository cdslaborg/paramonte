classdef SpecDRAM_ProposalStartCovMat_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_ProposalStartCovMat_class"
    end

    properties
        isPresent   = []
        Val         = []
        Def         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_ProposalStartCovMat_class(nd, methodName)
            self.isPresent  = false;
            self.Def        = eye(nd, nd);
            self.desc       = "ProposalStartCovMat is a real-valued positive-definite matrix of size (ndim,ndim), where ndim is the dimension of the "                  ...
                            + "sampling space. It serves as the best-guess starting covariance matrix of the proposal distribution. "                                   ...
                            + "To bring the sampling efficiency of " + methodName + " to within the desired requested range, the covariance matrix "                    ...
                            + "will be adaptively updated throughout the simulation, according to the user's requested schedule. If ProposalStartCovMat is not "        ...
                            + "provided by the user, its value will be automatically computed from the input variables ProposalStartCorMat and ProposalStartStdVec. "   ...
                            + "The default value of ProposalStartCovMat is an ndim-by-ndim Identity matrix."                                                            ...
                            ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, ProposalStartCovMat)
            if isempty(ProposalStartCovMat)
                self.Val        = self.Def;
            else
                self.Val        = ProposalStartCovMat;
                self.isPresent  = true;
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
                Err.msg         = Err.msg                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                        ...
                                + "The input requested proposalStartCovMat for the proposal of " + methodName   ...
                                + " is not a positive-definite matrix." + newline + newline                     ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end