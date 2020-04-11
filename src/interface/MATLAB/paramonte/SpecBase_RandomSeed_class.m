classdef SpecBase_RandomSeed_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_RandomSeed_class"
    end

    properties
        seed        = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

		function self = SpecBase_RandomSeed_class(methodName)
            self.seed 	= [];
            self.desc	= "randomSeed is a scalar 32bit integer that serves as the seed of the random number generator. When it is provided, "      ...
                        + "the seed of the random number generator will be set in a specific deterministic manner to enable future replications "   ...
                        + "of the simulation with the same configuration and input specifications. The default value for randomSeed is an integer " ...
                        + "vector of processor-dependent size and value that will vary from one simulation to another. "                            ...
                        + "However, enough care has been taken to assign unique random seed values to the random number generator on "              ...
                        + "each of the parallel threads (or images, processors, cores, ...) at all circumstances."                                  ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, randomSeed)
            if isempty(randomSeed)
                rng('shuffle');
            else
                self.seed = randomSeed;
                rng(abs(self.seed));
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end