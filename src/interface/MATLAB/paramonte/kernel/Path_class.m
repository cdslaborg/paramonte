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

classdef Path_class < handle

    properties (Constant)
        CLASS_NAME              = "@Path_class"
        MAX_FILE_PATH_LEN       = 2047
        WINDOWS_RESERVED_CHAR   = "<>:" + '"' + "|?*"
        SHELL_ESCAPE_CHAR       = " " ... — space character
                                + "!" ... — history expansion.
                                + '"' ... — shell syntax.
                                + "#" ... — comment start when preceded by whitespace; zsh wildcards.
                                + "$" ... — shell syntax.
                                + "&" ... — shell syntax.
                                + "'" ... — shell syntax.
                                + "(" ... — even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.
                                + ")" ... - even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.
                                + "*" ... — sh wildcard.
                                + "," ... — only inside brace expansion.
                                + ";" ... — shell syntax.
                                + "<" ... — shell syntax.
                                + "=" ... — in zsh, when it's at the beginning of a file name (filename expansion with PATH lookup).
                                + ">" ... — shell syntax.
                                + "?" ... — sh wildcard.
                                + "[" ... — sh wildcard.
                                + '\' ... — shell syntax.
                                + "]" ... — you may get away with leaving it unquoted.
                                + "^" ... — history expansion; zsh wildcard.
                                + "`" ... — shell syntax.
                                + "{" ... — brace expansion.
                                + "|" ... — shell syntax.
                                + "}" ... — needs to be escaped in zsh, other shells are more lenient when there's no matching opening brace.
                                + "~" ... — home directory expansion when at the beginning of a filename; zsh wildcard; safe when it's the last character.
                                ;
    end

    properties
        original                = []
        modified                = []
        dir                     = []
        name                    = []
        base                    = []
        ext                     = []
        slashOS                 = []
        Err                     = Err_class()
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = Path_class(inputPath, OS)
            self.queryPath(inputPath, OS);
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        % Construct a Path object as output. On output, check Err.occurred before using the output Path.

        function queryPath(self, inputPath, OS)

            FUNCTION_NAME = "@queryPath()";

            if ~isempty(inputPath)
                self.original       = strtrim(inputPath);
            else
                self.original       = "";
            end

            if ~isempty(OS)
                self.slashOS        = OS.slash;
                OSisWindows         = OS.isWindows;
            else
                OS = OS_class();
                OS.queryOS();
                if OS.Err.occurred
                    self.Err.occurred   = OS.Err.occurred;
                    self.Err.msg        = FUNCTION + ": Error occurred while querying OS type." + newline + self.Err.msg;
                end
                self.slashOS        = OS.slash;
                OSisWindows         = OS.isWindows;
                if self.Err.occurred, return; end
            end

            if OSisWindows
                self.modified       = self.winifyPath(self.original);
                if self.Err.occurred
                    self.Err.msg    = FUNCTION_NAME + ": Error occurred while making path='" + self.original + "' compatible with Windows OS." + newline + self.Err.msg;
                    return
                end
            else
                % if the path contains both / and \, then assume that it is already in linux style
                if index(self.original,"/") == 0   % path is given in Windows style
                    self.modified   = self.linifyPath(self.original);
                else
                    self.modified   = self.original;
                end
            end

            [self.dir, self.name, self.ext] = fileparts(self.modified);

            self.base = self.dir + self.slashOS + self.name;

        end % queryPath

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods(Static)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        % This code assumes that the input path is a Linux path. Windows paths like .\(paramonte)\paramonte.nml will be horribly
        % treated by this routine as \( also represents a Linux escape character. The result will be .(paramonte)\paramonte.nml
        % this routine strictly assumes that there is no dangling \ in the input Linux path, and if there is, then either it is used
        % to escape the special shell characters, or otherwise, the path is a Windows path.

        function outputPath = winifyPath(inputPath)

            FUNCTION_NAME   = "@winifyPath()";

            % note that multiple \ character in sequence is meaningless in Linux (basically \\ reduces to \),
            % and in Windows means the same as a single \. Therefore, reduce all sequential \ characters to a single \.
            outputPath = strtrim(inputPath);
            if ~isempty(inputPath)
                while true
                    outputPathNew   = strrep(outputPath, '\\', '\');
                    if outputPathNew == outputPath, break; end
                    outputPath      = outputPathNew;
                end
            end

            % Now check for the presence of any Linux Shell Escape Character in the input path without a preceding \.
            % If there is any, this would imply that the input path is a Windows path,
            % otherwise a escape character without preceding \ would be invalid in Linux.
            for i = 1 : length(Path_class.SHELL_ESCAPE_CHAR)
                if Path_class.SHELL_ESCAPE_CHAR(i:i) ~= '\'
                    for j = 1 : length(outputPath)
                        if outputPath(j:j) == Path_class.SHELL_ESCAPE_CHAR(i:i)
                            if j == 1
                                return  % it is a windows path, so no need for further winifying
                            elseif outputPath(j-1 : j-1) ~= '\'
                                return  % it is a windows path, so no need for further winifying
                            end
                        end
                    end
                end
            end

            % By now, there is no way but to assume that the path is indeed a Linux path
            % Thus, correct for any Linux Shell Escape Character in the input path:
            for i = 1 : length(Path_class.SHELL_ESCAPE_CHAR)
                outputPath = strrep(outputPath, '\' + Path_class.SHELL_ESCAPE_CHAR(i:i), Path_class.SHELL_ESCAPE_CHAR(i:i));
            end

            % check if the file name contains white space. if so, put the entire name in quotations
            function index_val = index(str, sub_str)
                index_val = strfind(str, sub_str);
                if isempty(index_val), index_val = 0; else, index_val = index_val(1); end
            end

            if index(outputPath, " ") ~= 0
                outputPath = '"' + outputPath  + '"';
            end
            outputPath = strrep(outputPath, '/', '\');

        end % winifyPath

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function outputPath = linifyPath(inputPath)

            FUNCTION_NAME   = "@linifyPath()";

            % check if the path is sandwiched between quotation marks. If so, remove them:

            outputPath      = strtrim(char(inputPath));
            i = length(outputPath);
            if i==0, return; end
            if i > 1
                if ((outputPath(1:1) == '"') && (outputPath(i:i)=='"')) || ((outputPath(1:1) == "'") && (outputPath(i:i) == "'"))
                    outputPathNew = outputPath(2 : i-1);
                else
                    outputPathNew = outputPath;
                end
            end

            % First change all backslashes to forward slash:
            outputPath = strrep(outputPathNew, '\', '/');

            % Now correct for any Linux Shell Escape Character in the input path:
            for i = 1 : length(Path_class.SHELL_ESCAPE_CHAR)
                if Path_class.SHELL_ESCAPE_CHAR(i:i) ~= '\'
                    outputPathNew   = strrep(outputPath, Path_class.SHELL_ESCAPE_CHAR(i:i), '\' + Path_class.SHELL_ESCAPE_CHAR(i:i));
                    outputPath      = outputPathNew;
                end
            end

        end % linifyPath

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [Err, slashOS] = getSlashOS(slashOS)
            FUNCTION_NAME   = "@getSlashOS()";
            Err             = Err_class();
            OS              = OS_class();
            OS.queryOS();

            if OS.Err.occurred
                Err     = OS.Err;
                Err.msg = FUNCTION_NAME + ": Error occurred while fetching the OS slash character." + newline + Err.msg;
                return
            end

            if OS.isWindows
                slashOS = '\';
            else
                slashOS = '/';
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [Err, outputPath] = modifyPath(inputPath)
            FUNCTION_NAME   = "@modifyPath()";
            Err             = Err_class();
            outputPath      = strtrim(inputPath);
            OS              = OS_class();
            OS.queryOS();

            if OS.Err.occurred
                Err         = OS.Err;
                Err.msg     = FUNCTION_NAME + ": Error occurred while modifying inputPath='" + outputPath + "'." + newline + Err.msg;
                return
            end

            if OS.isWindows
                outputPath  = Path_class.winifyPath(inputPath);
                if Err.occurred
                    Err.msg =  FUNCTION_NAME + ": Error occurred while making path='"       ...
                            + inputPath + "' compatible with Windows OS." + newline + Err.msg;
                    return
                end
            else
                outputPath  = Path_class.linifyPath(inputPath);
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end