classdef SpecMCMC_ChainSize_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecMCMC_ChainSize_class"
    end

    properties
        val         = []
        def         = []
        null        = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecMCMC_ChainSize_class(methodName)
            self.def    = 100000;
            self.desc   = "chainSize determines the number of non-refined, potentially auto-correlated, but unique, samples "                   ...
                        + "drawn by the MCMC sampler before stopping " + methodName + ". "                                                      ...
                        + "For example, if you specify chainSize = 10000, then 10000 unique sample points (with no duplicates) will be drawn "  ...
                        + "from the target objective function that the user has provided. "                                                     ...
                        + "The input value for chainSize must be a positive integer of a minimum value ndim+1 or larger, "                      ...
                        + "where ndim is the number of variables that define the domain of the objective function to be sampled. "              ...
                        + "The default value is " + num2str(self.def) + "."                                                                     ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, chainSize)
            self.val = chainSize;
            if isempty(self.val)
                self.val = self.def;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName, nd)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < nd + 1
                Err.occurred    = true;
                Err.msg         =   Err.msg                                                                                         ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                            ...
                                + "The input requested value for chainSize (" + num2str(self.val) + ") can neither be negative "    ...
                                + "nor smaller than ndim+1, where ndim represents the dimension of the sampling space, here ndim="  ...
                                + num2str(nd) + ". "                                                                                ...
                                + "If you don't know an appropriate value for chainSize, drop it from the input list. "             ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline         ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
