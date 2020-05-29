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