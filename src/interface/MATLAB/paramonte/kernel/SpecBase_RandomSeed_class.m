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