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
