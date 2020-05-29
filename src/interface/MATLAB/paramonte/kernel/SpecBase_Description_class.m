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
classdef SpecBase_Description_class < handle

    properties (Constant)
        CLASS_NAME          = "@SpecBase_Description_class"
        MAX_DESCRIPTION_LEN = 4096
    end

    properties
        val                 = []
        def                 = []
        desc                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_Description_class(methodName)
            self.def    = "Nothing provided by the user.";
            self.desc   = "The variable 'description' contains general information about the specific " + methodName + " simulation that "              ...
                        + "is going to be performed. It has no effects on the simulation and serves only as a general description of the simulation "   ...
                        + "for future reference. The " + methodName + " parser automatically recognizes the C-style '\\n' escape sequence as the "      ...
                        + "new-line character, and '\\\\' as the backslash character '\\' if they used in the description. For example, '\\\\n' "       ...
                        + "will be converted to '\\n' on the output, while '\\n' translates to the new-line character. Other C escape sequences "       ...
                        + "are neither supported nor needed. The default value for description is '" + self.def + "'."                                  ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, description)
        
            if isempty(description)
                self.val = strtrim(self.def);
            else
                self.val = strtrim(description);
            end
            
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%**********************************************************************************************************************************

end