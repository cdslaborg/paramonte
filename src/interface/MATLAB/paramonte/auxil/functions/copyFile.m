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
function Err = copyFile(pathOld, pathNew, isWindows)
    
    FUNCTION_NAME   = "@copyFile()";
    
    Err             = Err_class();
    Err.occurred    = false;

    if length(strtrim(pathOld)) == 0, return; end

    % First check whether file exists:
    
    fileExists      = isfile(pathNew);   % check if the file already exists
    
    if fileExists
        Err.occurred    = true;
        Err.msg         = FUNCTION_NAME + ": The requested copy file = '" + pathNew + "' already exists.";
        return
    end

    if isWindows
        cmd = 'copy "'  + pathOld + '" "'   + pathNew + '" > nul';
    else
        cmd = "cp "     + pathOld + " "     + pathNew;
    end
    
    counter = 0;
    while true
        counter = counter + 1;
        Err     = SystemInfo_class.executeCmd(cmd);
        if Err.occurred
            Err.msg     = FUNCTION_NAME + ": Error occurred while executing command " + cmd + "'.\n";
            return
        end
        % ensure file is copied
        fileExists      = isfile(pathNew);   % check if the file already exists
        if ~fileExists && counter < 100, continue; end
        break
    end
    
    if ~fileExists
        Err.occurred   = true;
        Err.msg        = FUNCTION_NAME + ": Failed to copy file from '" + pathOld + "' to '" + pathNew + "' after " + num2str(counter) + " attempts.";
        return
    end
    
end