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
classdef SpecBase_OutputFileName_class < Path_class

    properties (Constant)
        SUB_CLASS_NAME = "@SpecBase_OutputFileName_class"
    end

    properties
        def         = []
        namePrefix  = []
        pathPrefix  = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_OutputFileName_class(methodName)
            
            self = self@Path_class([], []);
            
            self.def    = methodName + '_run_' + datestr(now, 'ddmmyy') + '_' + datestr(now, 'HHMMSS_FFF');
            self.desc   = "outputFileName contains the path and the base of the filename for " + methodName + " output files. "                         ...
                        + "If not provided by the user, the default outputFileName is constructed from the current date and time:"                      ...
                        + newline + newline + Decoration_class.TAB + methodName + "_run_ddmmyy_hhmmss_fff" + newline + newline                          ...
                        + "where dd, mm, yy, hh, mm, ss stand respectively for the current day, month, year, hour, minute, second. "                    ...
                        + "In such a case, the default directory for the output files will be the current working directory of "                        ...
                        + methodName + ". If outputFileName is provided, but ends with a separator character '/' or '\' (as in Linux or Windows OS), "  ...
                        + "then its value will be used as the directory to which " + methodName + " output files will be written. In this case, "       ...
                        + "the output file naming convention described above will be used. Also, the given directory will be automatically created "    ...
                        + "if it does not exist already."                                                                                               ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, outputFileName)
            if isempty(outputFileName)
                self.original = self.def;
            else
                self.original = strtrim(outputFileName);
            end
            self.modified = self.original;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end