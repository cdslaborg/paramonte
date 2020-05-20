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
function Err = removeFile(path,isWindows)

    FUNCTION_NAME   = "@removeFile()";

    Err             = Err_class();
    Err.occurred    = false;

    % First check whether file exists:

    fileExists      = isfile(path); % check if the file already exists

    if ~fileExists
        Err.occurred    = true;
        Err.msg         = FUNCTION_NAME + ": The requested file = '" + path + "' does not exist.";
        return
    end

    if isWindows
        cmd = "del " + path + " > nul";
    else
        cmd = "rm " + path;
    end

    counter = 0;
    while true
        counter = counter + 1;
        Err     = executeCmd(cmd);
        if Err.occurred
            Err.msg     = FUNCTION_NAME + ": Error occurred while executing command " + cmd + "'.\n";
            return
        end

        % ensure file is removed

        fileExists  = isfile(path); % check if the file already exists
        if fileExists
            Err.occurred    = true;
            Err.msg         = FUNCTION_NAME + ": Failed to remove file = '" + path + "' after " + num2str(counter) + " attempts.";
            return
        end
        if fileExists && counter < 100, continue; end
        break
    end

    if fileExists
        Err.occurred    = true;
        Errmsg          = FUNCTION_NAME + ": Failed to remove file = '" + path + "' after " + num2str(counter) + " attempts.";
        return
    end

end