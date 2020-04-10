classdef SpecDRAM_AdaptiveUpdateCount_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_AdaptiveUpdateCount_class"
    end

    properties
        val         = []
        def         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_AdaptiveUpdateCount_class(methodName)
            self.def    = Constants.POSINF_IK;
            self.desc   = "adaptiveUpdateCount represents the total number of adaptive updates that will be made "                              ...
                        + "to the parameters of the proposal distribution, to increase the efficiency of the sampler "                          ...
                        + "thus increasing the sampling efficiency of " + methodName + ". "                                                     ...
                        + "Every adaptiveUpdatePeriod number of calls to the objective function, the parameters of the proposal distribution "  ...
                        + "will be updated until either the total number of adaptive updates reaches the value of adaptiveUpdateCount. "        ...
                        + "This variable must be a non-negative integer. As a rule of thumb, it may be appropriate to "                         ...
                        + "set the input variable chainSize > 2 * adaptiveUpdatePeriod * adaptiveUpdateCount, "                                 ...
                        + "to ensure ergodicity and stationarity of the MCMC sampler. "                                                         ...
                        + "If adaptiveUpdateCount=0, then the proposal distribution parameters will be fixed to "                               ...
                        + "the initial input values throughout the entire MCMC sampling. "                                                      ...
                        + "The default value is " + num2str(self.def) + "."                                                                     ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, adaptiveUpdateCount)
            if isempty(adaptiveUpdateCount)
                self.val = self.def;
            else
                self.val = adaptiveUpdateCount;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < 0
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                               ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                                ...
                                + "The input requested value for adaptiveUpdateCount (" + num2str(self.val) + ") "                                      ...
                                + "can not be negative. If you are not sure of the appropriate value for adaptiveUpdateCount, drop it "                 ...
                                + "from the input list. " + methodName + " will automatically assign an appropriate value to it." + newline + newline   ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end