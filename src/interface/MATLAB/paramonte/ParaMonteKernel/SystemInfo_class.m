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

            Err = executeCmd("del " + filename);
            if (Err.occurred == 0) && (Err.msg == ("Could Not Find " + pwd + '\' + filename + newline))
                Err.msg =   FUNCTION_NAME + ": Error occurred while attempting to delete file " + filename;
                Err.abort();
                return
            end

            Err = executeCmd("del " + stdErr);
            if (Err.occurred == 0) && (Err.msg == ("Could Not Find " + pwd + '\' + stdErr + newline))
                Err.msg =   FUNCTION_NAME + ": Error occurred while attempting to delete file " + stdErr;
                Err.abort();
                return
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end