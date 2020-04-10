classdef SpecDRAM_DelayedRejectionScaleFactorVec_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_DelayedRejectionScaleFactorVec_class"
    end

    properties
        Val         = []
        def         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_DelayedRejectionScaleFactorVec_class(nd, methodName)
            self.def    = 0.5^(1/nd);  % This gives a half volume to the covariance
            self.desc   = "delayedRejectionScaleFactorVec is a real-valued positive vector of length (1:delayedRejectionCount) by which "                       ...
                        + "the covariance matrix of the proposal distribution of " + methodName + " sampler is scaled when the Delayed Rejection (DR) "         ...
                        + "scheme is activated (by setting delayedRejectionCount>0). At each ith stage of the DR process, "                                     ...
                        + "the proposal distribution from the last stage is scaled by the factor delayedRejectionScaleFactorVec(i). "                           ...
                        + "Missing elements of the delayedRejectionScaleFactorVec in the input to " + methodName + " will be set to the default value. "        ...
                        + "The default value at all stages is 0.5^(1/ndim) = " + num2str(self.def) + ", which reduces the "                                     ...
                        + "volume of the covariance matrix of the proposal from the last DR stage by one half. The variable ndim represents the "               ...
                        + "number of dimensions of the Domain of the objective function."                                                                       ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, delayedRejectionScaleFactorVec, delayedRejectionCount)
            if ~isempty(delayedRejectionScaleFactorVec)
                for i = 1 : length(delayedRejectionScaleFactorVec)
                    self.Val(i) = delayedRejectionScaleFactorVec(i);
                end
            end
            if isempty(self.Val) && (delayedRejectionCount > 0)
                for i = 1 : delayedRejectionCount
                    self.Val(i) = self.def;
                end
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, delayedRejectionCount)
            FUNCTION_NAME = "@checkForSanity()";
            if length(self.Val) ~= delayedRejectionCount
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                                        ...
                                + "The length of the vector delayedRejectionScaleFactorVec (" + num2str(length(self.Val)) + ")"                                 ...
                                + " is not equal to delayedRejectionCount = " + num2str(delayedRejectionCount)                                                  ...
                                + ". If you are not sure how to set the values of delayedRejectionScaleFactorVec, drop it from the input. "                     ...
                                + methodName + " will automatically set the appropriate value for delayedRejectionScaleFactorVec." + newline + newline          ...
                                ;
            end
            for i = 1 : length(self.Val)
                if self.Val(i) <= 0
                    Err.occurred    = true;
                    Err.msg         = Err.msg                                                                                                                   ...
                                    + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. " + "The input value for the element "                               ...
                                    + num2str(i) + " of the variable delayedRejectionScaleFactorVec cannot be smaller than or equal to 0." + newline + newline  ...
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