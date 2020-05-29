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
classdef SpecDRAM_ProposalModel_class < handle

    properties (Constant)
        CLASS_NAME              = "@SpecDRAM_ProposalModel_class"
        MAX_LEN_PROPOSAL_MODEL  = 63
    end

    properties
        isUniform               = []
        isNormal                = []
        uniform                 = []
        normal                  = []
        val                     = []
        def                     = []
        desc                    = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_ProposalModel_class()
            self.isUniform  = false;
            self.isNormal   = false;
            self.uniform    = "uniform";
            self.normal     = "normal";
            self.def        = self.normal;
            self.desc       = "proposalModel is a string variable containing the name of the proposal distribution for the MCMC sampler. "              ...
                            + "The string value must be enclosed by either single or double quotation marks when provided as input. "                   ...
                            + "One option is currently supported:" + newline + newline                                                                  ...
                            + "    proposalModel = '" + self.uniform + "'" + newline + newline                                                          ...
                            + "            The proposals will be drawn uniformly from within a ndim-dimensional ellipsoid whose covariance matrix "     ...
                            + "and scale are initialized by the user and optionally adaptively updated throughout the simulation." + newline + newline  ...
                            + "The default value is '" + self.def + "'."                                                                                ...
                            ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, proposalModel)
            if isempty(proposalModel)
                self.val = self.def;
            else
                self.val = lower(strtrim(proposalModel));
            end
            if self.val == self.normal,     self.isNormal   = true; end
            if self.val == self.uniform,    self.isUniform  = true; end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if ~(self.isNormal || self.isUniform)
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                        ...
                                + "Invalid requested value for the proposalModel of " + methodName + ". The input requested "   ...
                                + "proposal model (" + self.val + ") is not supported. "                                        ...
                                + "The variable proposalModel cannot be set to anything other than '"                           ...
                                + self.normal + "', or '" + self.uniform + "'." + newline + newline                             ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end