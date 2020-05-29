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
classdef SpecBase_TargetAcceptanceRate_class < handle

    properties (Constant)
        CLASS_NAME          = "@SpecBase_TargetAcceptanceRate_class"
    end

    properties
        scalingRequested    = []
        val                 = []
        desc                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        function self = SpecBase_TargetAcceptanceRate_class(methodName)
            self.scalingRequested   = true;
            self.desc               = "targetAcceptanceRate sets an optimal target for the ratio of the number of accepted objective function calls to the "    ...
                                    + "total number of function calls by " + methodName + ". By default, it is a real number between 0 and 1. "                 ...
                                    + "If provided by the user, " + methodName + " will attempt (but not guarantee) to bring the average acceptance ratio "     ...
                                    + "of the sampler as close to the user-provided target ratio as possible. The success of " + methodName                     ...
                                    + " in keeping the average acceptance ratio close to the requested target value depends heavily on:" + newline              ...
                                    + "    1) the value of adaptiveUpdatePeriod; the larger, the easier." + newline                                             ...
                                    + "    2) the value of adaptiveUpdateCount; the larger, the easier." + newline                                              ...
                                    + "Note that the acceptance ratio adjustments will only occur every adaptiveUpdatePeriod sampling steps for a total "       ...
                                    + "number of adaptiveUpdateCount. "                                                                                         ...
                                    + "There is no default value for targetAcceptanceRate, as the acceptance ratio is not directly adjusted during sampling."   ...
                                    ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, targetAcceptanceRate)
            if isempty(targetAcceptanceRate)
                self.scalingRequested = false;
            else
                self.val = targetAcceptanceRate;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err)
            FUNCTION_NAME = "@checkForSanity()";
            if ~self.scalingRequested, return; end
            if self.val <= 0
                Err.occurred    = true;
                Err.msg         = Err.msg...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                    ...
                                + "The target acceptance ratio targetAcceptanceRate (" + num2str(self.val)  ...
                                + ") cannot be less than or equal to 0." + newline + newline                ...
                                ;
            end
            if self.val >= 1
                Err.occurred    = true;
                Err.msg         = Err.msg...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                    ...
                                + "The target acceptance ratio targetAcceptanceRate (" + num2str(self.val)  ...
                                + ") cannot be larger than or equal to 1." + newline + newline              ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end