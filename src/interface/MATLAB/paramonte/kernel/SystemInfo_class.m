%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        function self = SystemInfo_class(isWindowsOS, systemInfoFilePath)
            
            if nargin==2

                % first check and see if any cache file containing the system info exists or not

                self.Err     = Err_class();
                self.Err.msg = "";
                self.Err.occurred = false;
                if isfile(systemInfoFilePath)
                    self.info = fileread(systemInfoFilePath);
                    self.info = strrep(self.info, '\', '\\');
                    return
                else
                    self.Err.occurred = true;
                end

            else

                [self.info, self.Err, self.nRecord] = self.getSystemInfo(isWindowsOS);

            end

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