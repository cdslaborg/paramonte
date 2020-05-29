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