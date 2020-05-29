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
classdef SpecDRAM_AdaptiveUpdatePeriod_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecDRAM_AdaptiveUpdatePeriod_class"
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

        function self = SpecDRAM_AdaptiveUpdatePeriod_class(nd,methodName)
            self.def    = nd * 4;   % max(nd+1,100)
            self.desc   = "Every adaptiveUpdatePeriod calls to the objective function, the parameters of the proposal distribution will be updated. "   ...
                        + "In parallel mode, this update will happen every adaptiveUpdatePeriod or more calls to the objective function. "              ...
                        + "The variable adaptiveUpdatePeriod must be a positive integer of value ndim+1 or larger. "                                    ...
                        + "The larger the value of adaptiveUpdatePeriod, the easier it will be "                                                        ...
                        + "for " + methodName + " kernel to keep the sampling efficiency within the desired range (if requested). "                     ...
                        + "However, too large values for adaptiveUpdatePeriod will only delay the adaptation of the proposal distribution to "          ...
                        + "the global structure of the objective function that is being sampled. "                                                      ...
                        + "If adaptiveUpdatePeriod>=chainSize, then no adaptive updates to the proposal distribution will be made. "                    ...
                        + "The default value is 4 * ndim, where ndim is the dimension of the sampling space. "                                   ...
                        + "In this particular " + methodName + " simulation, this corresponds to the value "                                            ...
                        + num2str(self.def) + "."                                                                                                       ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, adaptiveUpdatePeriod)
            if isempty(adaptiveUpdatePeriod)
                self.val = self.def;
            else
                self.val = adaptiveUpdatePeriod;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < 1
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                        ...
                                + "Invalid requested value for adaptiveUpdatePeriod. "                                                          ...
                                + "The input requested value for adaptiveUpdatePeriod (" + num2str(self.val)                                    ...
                                + ") cannot be less than 1. If you are not sure of the appropriate value for adaptiveUpdatePeriod, drop it "    ...
                                + "from the input list. " + methodName + " will automatically assign an appropriate value to it."               ...
                                + newline + newline                                                                                             ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end