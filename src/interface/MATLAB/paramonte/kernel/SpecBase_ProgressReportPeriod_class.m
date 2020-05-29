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
classdef SpecBase_ProgressReportPeriod_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_ProgressReportPeriod_class"
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

        function self = SpecBase_ProgressReportPeriod_class()
            self.def 	= 1000;
            self.desc	= "Every progressReportPeriod calls to the objective function, the sampling progress will be reported to the log file. "...
						+ "Note that progressReportPeriod must be a positive integer. "                                                         ...
						+ "The default value is " + num2str(self.def) + "."                                                                     ...
						;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, progressReportPeriod)
            self.val = progressReportPeriod;
            if isempty(self.val)
                self.val = self.def;
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
                                + "The input value for variable progressReportPeriod must be a positive integer value. If you are not sure "    ...
                                + "about the appropriate value for this variable, simply drop it from the input. "                              ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline                     ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end