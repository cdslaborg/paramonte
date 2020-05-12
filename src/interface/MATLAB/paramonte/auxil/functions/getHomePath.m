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
%   we ask you to acknowledge the ParaMonte library's usage
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function homePath = getHomePath(varargin)
    freshRequested = false;
    persistent homePathPersistent
    if nargin==0
        if isempty(homePathPersistent); freshRequested = true; end
    elseif nargin==1 && ( strcmpi(string(varargin{1}),"fresh") || strcmpi(string(varargin{1}),"new") || strcmpi(string(varargin{1}),"reset") )
        freshRequested = true;
    else
        error( "getHomePath() takes at most one argument of the following values: ""new"", ""fresh"", ""reset"", all with the same meaning." );
    end
    if freshRequested
        if ispc
            [errStat,homePathPersistent] = system("echo %HOMEPATH%");
            if errStat  
                error   ( newline ...
                        + "Failed to capture the path to the home directory of the Windows system. This is highlt unusual and likely a low-level MATLAB problem." ...
                        + newline ...
                        );
            else
                homePathPersistent = strtrim(homePathPersistent);
            end
        else
            homePathPersistent = strtrim(getFullPath("~"));
        end
    end
    homePath = string(homePathPersistent);
end