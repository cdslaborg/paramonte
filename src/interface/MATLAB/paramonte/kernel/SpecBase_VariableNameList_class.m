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
classdef SpecBase_VariableNameList_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_VariableNameList_class"
    end

    properties
        Val         = []
        Def         = []
        desc        = []
        prefix      = []
        MaxLen      = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_VariableNameList_class(methodName, nd)
            self.prefix = "SampleVariable";
            self.Val    = cell(nd,1);
            self.Def    = cell(nd,1);
            for i = 1 : nd
                self.Def{i} = strtrim(self.prefix + num2str(i));
            end
            self.desc   = "variableNameList contains the names of the variables to be sampled by " + methodName + ". "...
                        + "It is used to construct the header of the output sample file. "...
                        + "Any element of variableNameList that is not set by the user will be automatically assigned a default name. "...
                        + "The default value is '" + self.prefix + "i' where integer 'i' is the index of the variable."...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, variableNameList)
            self.MaxLen.val = -1;
            if isempty(variableNameList)
                self.Val = self.Def;
                self.MaxLen.val = 15;
            else
                if length(variableNameList) <= length(self.Val)
                    for i = 1 : length(variableNameList)
                        self.Val{i} = variableNameList(i);
                        lentrim = length(char(strtrim(self.Val{i})));
                        if lentrim > self.MaxLen.val, self.MaxLen.val = lentrim; end
                    end
                    for i = length(variableNameList)+1 : length(self.Val)
                        self.Val{i} = self.Def(i);
                    end
                else
                    for i = 1 : length(self.Val)
                        self.Val{i} = variableNameList(i);
                        lentrim = length(char(strtrim(self.Val{i})));
                        if lentrim > self.MaxLen.val, self.MaxLen.val = lentrim; end
                    end
                end
            end
            self.MaxLen.str = num2str(self.MaxLen.val);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end