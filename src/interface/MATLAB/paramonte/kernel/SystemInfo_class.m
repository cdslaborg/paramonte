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
classdef SystemInfo_class < handle

    properties (Constant)
        CLASS_NAME      = "@System_class"
        MAX_OS_NAME_LEN = 63
    end

    properties
        nRecord         = []
        List            = []
        Err             = []
        info            = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SystemInfo_class(isWindowsOS)
            [self.info, self.Err, self.nRecord] = self.getSystemInfo(isWindowsOS);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [info, Err, count] = getSystemInfo(isWindowsOS)

            FUNCTION_NAME = SystemInfo_class.CLASS_NAME + "@getSystemInfo()";
            Err     = Err_class();
            Err.msg = "";

            % generate a brand new, non-existing filename

            RFN = RandomFileName_class([], "getFileList", []);
            if RFN.Err.occurred
                RFN.Err.msg = FUNCTION_NAME + RFN.Err.msg;
                return
            end
            filename = RFN.path;
            stdErr = filename + ".stderr";

            % determine the operating system

            if ~isempty(isWindowsOS)
                isWindows = isWindowsOS;
            else
                OS = OS_class();
                OS.queryOS();
                if OS.Err.occurred
                    Err     = OS.Err;
                    Err.msg = FUNCTION_NAME + Err.msg;
                    return
                end
                isWindows = OS.isWindows;
            end

            if isWindows    % it's Windows cmd
                command = "systeminfo > " + filename;
            elseif ismac
                command = "uname -a >> " + filename + "; sysctl -a | grep machdep.cpu >> " + filename;
            else
                command = "uname -a >> " + filename + "; lscpu >> " + filename;
            end

            Err = executeCmd(command + " 2> " + stdErr);
            if Err.occurred
                Err.msg =   FUNCTION_NAME + ": Error occurred while attempting to write the system info to external file." + newline + Err.msg;
                disp(Err.msg);
                return
            end

            % now read the contents of the file

            info = fileread(filename);
            info = strrep(info, '\', '\\');

            count = [];

            % remove the files

            if exist(filename)
                delete(filename);
            else
%                Err.msg = FUNCTION_NAME + "Could Not Find " + pwd + '/' + filename + newline;
%                Err.abort();
%                return
            end

            if exist(stdErr)
                delete(stdErr);
            else
%                Err.msg = FUNCTION_NAME + "Could Not Find " + pwd + '/' + filename + newline;
%                Err.abort();
%                return
            end

%            Err = executeCmd("del " + filename);
%            if (Err.occurred == 0) && (Err.msg == ("Could Not Find " + pwd + '\' + filename + newline))
%                Err.msg =   FUNCTION_NAME + ": Error occurred while attempting to delete file " + filename;
%                Err.abort();
%                return
%            end
%
%            Err = executeCmd("del " + stdErr);
%            if (Err.occurred == 0) && (Err.msg == ("Could Not Find " + pwd + '\' + stdErr + newline))
%                Err.msg =   FUNCTION_NAME + ": Error occurred while attempting to delete file " + stdErr;
%                Err.abort();
%                return
%            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end