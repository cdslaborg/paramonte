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
function out = parseArgs(self,varargin)

    vararginLen = length(varargin);
    if isa(self,"struct")
        selfProperties = fieldnames(self);
        selfPropertiesLen = length(selfProperties);
    else
        selfProperties = properties(self);
        selfPropertiesLen = length(selfProperties);
    end

    for i = 1:2:vararginLen
        propertyDoesNotExist = true;
        vararginString = string(varargin{i});
        for ip = 1:selfPropertiesLen
            if strcmp(vararginString,string(selfProperties(ip)))
                propertyDoesNotExist = false;
                if i < vararginLen
                    if isa(self.(selfProperties{ip}),"struct") && isa(varargin{i+1},"cell")
                        self.(selfProperties{ip}) = parseArgs( self.(selfProperties{ip}) , varargin{i+1}{:} );
                    else
                        self.(selfProperties{ip}) = varargin{i+1};
                    end
                else
                    error("The corresponding value for the property """ + string(selfProperties{ip}) + """ is missing as input argument.");
                end
                break;
            end
        end
        if propertyDoesNotExist
            error("The requested property """ + string(varargin{i}) + """ does not exist.");
        end
    end

    out = self;

end
