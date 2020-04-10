classdef SpecDRAM_GreedyAdaptationCount_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_GreedyAdaptationCount_class"
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

        function self = SpecDRAM_GreedyAdaptationCount_class(methodName)
            self.def    = 0;
            self.desc   = "If greedyAdaptationCount is set to a positive integer then the first greedyAdaptationCount number of "                       ...
                        + "the adaptive updates of the sampler will be made using only the 'unique' accepted points in the MCMC chain. "                ...
                        + "This is useful for example, the function to be sampled by " + methodName + " is high dimensional, "                          ...
                        + "in which case, the adaptive updates to " + methodName + "'s sampler distribution will less likely lead to "                  ...
                        + "numerical instabilities, for example, a singular covariance matrix for the multivariate proposal sampler. "                  ...
                        + "The variable greedyAdaptationCount must be a non-negative integer, and not larger than the value of adaptiveUpdateCount. "   ...
                        + "If it is larger, it will be automatically set to adaptiveUpdateCount for the simulation. "                                   ...
                        + "The default value is " + num2str(self.def) + "."                                                                             ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, greedyAdaptationCount)
            if isempty(greedyAdaptationCount)
                self.val = self.def;
            else
                self.val = greedyAdaptationCount;
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
                                + "The input requested value for greedyAdaptationCount (" + num2str(self.val) + ") "                                    ...
                                + "can not be negative. If you are not sure of the appropriate value for greedyAdaptationCount, drop it "               ...
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