%>  \brief
%>  Return the Get absolute canonical path of a file or folder.
%>  Absolute path names are safer than relative paths, e.g. when
%>  a GUI or TIMER callback changes the current directory.
%>  Only canonical paths without ``"."`` and ``".."``
%>  can be recognized uniquely.<br>
%>  Long path names (>259 characters) require a magic initial
%>  key ``"\\?\"`` to be handled by Windows API functions, e.g. for
%>  MATLAB intrinsic routines ``fopen()``, ``dir()`` and ``exist()``.
%>
%>  \note
%>  Some functions of the Windows-API still do not support long file names.
%>  E.g. the Recycler and the Windows Explorer fail even with the magic ``'\\?\'`` prefix.
%>  Some functions of MATLAB accept 260 characters (value of MAX_PATH), some at 259 already.
%>  The ``'fat'`` style is useful e.g., when MATLAB ``dir()`` command is called for a folder
%>  with less than ``260`` characters, but together with the file name this limit is exceeded.
%>  Then, ``"dir(pm.sys.path.abs([folder, '\*.*'], 'fat'))"`` helps.
%>
%>  \devnote
%>  A MEX version of this function developed by Jan Simon performs much faster on Windows.
%>  Difference between M- and Mex-version:<br>
%>      - Mex does not work under MacOS/Unix.
%>      - Mex calls Windows API function pm.sys.path.abs.
%>      - Mex is much faster.
%>
%>  \param[in]  path    :   The input argument that can be either,<br>
%>                              -   a scalar MATLAB string or character, or,
%>                              -   a MATLAB cell array of strings or character values,
%>                          containing the absolute or relative name(s) of a file(s) or folder(s).
%>                          The path need not exist. Unicode strings, UNC paths and long
%>                          names are supported.
%>  
%>  \param[in]  style   :   The optional input scalar MATLAB string or character,
%>                          containing the style of the output as string.
%>                          The following values are possible:<br>
%>                              -   'auto'  :   Add '\\?\' or '\\?\UNC\' for long names on demand.
%>                              -   'lean'  :   Magic string is not added.
%>                              -   'fat'   :   Magic string is added for short names also.
%>                          The input ``style`` is ignored when not running under Windows.
%>                          (**optional**, default = 'auto')
%>
%>  \return
%>  `path`              :   The output scalar MATLAB character string
%>                          containing the absolute path corresponding to the input path.
%>                          If the input ``path`` is empty, ``pathAbs`` is the current directory.
%>                          The optional prefixes ``'\\?\'`` or ``'\\?\UNC'`` are added on demand
%>                          as requested by the optional input argument ``style``.
%>
%>  \interface{abs}
%>  \code{.m}
%>
%>      path = pm.sys.path.abs(file)
%>      path = pm.sys.path.abs(file, style)
%>
%>  \endcode
%>
%>  \example{abs}
%>
%>      cd(tempdir);                                        % Assumed as 'C:\Temp' here
%>      pm.sys.path.abs('File.Ext')                         % 'C:\Temp\File.Ext'
%>      pm.sys.path.abs('..\File.Ext')                      % 'C:\File.Ext'
%>      pm.sys.path.abs('..\..\File.Ext')                   % 'C:\File.Ext'
%>      pm.sys.path.abs('.\File.Ext')                       % 'C:\Temp\File.Ext'
%>      pm.sys.path.abs('*.txt')                            % 'C:\Temp\*.txt'
%>      pm.sys.path.abs('..')                               % 'C:\'
%>      pm.sys.path.abs('..\..\..')                         % 'C:\'
%>      pm.sys.path.abs('Folder\')                          % 'C:\Temp\Folder\'
%>      pm.sys.path.abs('D:\A\..\B')                        % 'D:\B'
%>      pm.sys.path.abs('\\Server\Folder\Sub\..\File.ext')  % '\\Server\Folder\File.ext'
%>      pm.sys.path.abs({'..', 'new'})                      % {'C:\', 'C:\Temp\new'}
%>      pm.sys.path.abs('.', 'fat')                         % '\\?\C:\Temp\File.Ext'
%>  See also: http://www.n-simon.de/mex - A collection of MATLAB mex files for system path manipulation.
%>
%>  \final{abs}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:15 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>      This file is based on a functionality originally developed by Jan Simon.
%>
%>      Copyright (c) 2016, Jan Simon
%>      All rights reserved.
%>
%>      Redistribution and use in source and binary forms, with or without
%>      modification, are permitted provided that the following conditions are met:
%>
%>      * Redistributions of source code must retain the above copyright notice, this
%>        list of conditions and the following disclaimer.
%>
%>      * Redistributions in binary form must reproduce the above copyright notice,
%>        this list of conditions and the following disclaimer in the documentation
%>        and/or other materials provided with the distribution
%>      * Neither the name of  nor the names of its
%>        contributors may be used to endorse or promote products derived from this
%>        software without specific prior written permission.
%>      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%>      AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%>      IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%>      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%>      FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%>      DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%>      SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%>      CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%>      OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%>      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%>
%>      Tested: MATLAB 2009a, 2015b(32/64), 2016b, 2018b, Win7/10
%>              Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
%>      Assumed Compatibility: higher MATLAB versions
%>      Author: Jan Simon, Heidelberg, (C) 2009-2019 matlab.2010(a)n(MINUS)simon.de
%>
%>      $JRev: R-M V:038 Sum:C/6JMzUYsYsc Date:19-May-2019 17:25:55 $
%>      $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
%>      $UnitTest: uTest_getFullPath $
%>      $File: Tools\GLFile\pm.sys.path.abs.m $
%>      History:
%>      001: 20-Apr-2010 22:28, Successor of Rel2AbsPath.
%>      010: 27-Jul-2008 21:59, Consider leading separator in M-version also.
%>      011: 24-Jan-2011 12:11, Cell strings, '~File' under linux.
%>           Check of input types in the M-version.
%>      015: 31-Mar-2011 10:48, BUGFIX: Accept [] as input as in the Mex version.
%>           Thanks to Jiro Doke, who found this bug by running the test function for
%>           the M-version.
%>      020: 18-Oct-2011 00:57, BUGFIX: Linux version created bad results.
%>           Thanks to Daniel.
%>      024: 10-Dec-2011 14:00, Care for long names under Windows in M-version.
%>           Improved the unittest function for Linux. Thanks to Paul Sexton.
%>      025: 09-Aug-2012 14:00, In MEX: Paths starting with "\\" can be non-UNC.
%>           The former version treated "\\?\C:\<longpath>\file" as UNC path and
%>           replied "\\?\UNC\?\C:\<longpath>\file".
%>      032: 12-Jan-2013 21:16, 'auto', 'lean' and 'fat' style.
%>      038: 19-May-2019 17:25, BUGFIX, Thanks HHang Li, "File(7:..." -> "File(8:..."
function path = abs(path, style)

    % Set the magic prefix for long Windows names.

    if nargin < 2
        style = 'auto';
    elseif isstring(style)
        style = convertStringsToChars(style);
    else
        error   ( newline ...
                + "The input argument `style` must be a string. See the function documentation for possible string values." + newline ...
                + newline ...
                );
    end

    if isstring(path)
        path = convertStringsToChars(path);
    elseif isa(path, 'cell')
        for icell = 1 : numel(path)
            path{icell} = pm.sys.path.abs(path{icell}, style);
        end
        return;
    end

    % Check this once only.

    isWindows = strncmpi(computer, 'PC', 2);
    MAX_PATH = 260;

    % Warn once per session (disable this under Linux/MacOS).

    persistent hasDataRead
    if isempty(hasDataRead)
        % Test this once only - there is no relation to the existence of DATAREAD!
        %if isWindows
        %   Show a warning, if the slower MATLAB version is used - commented, because
        %   this is not a problem and it might be even useful when the MEX-folder is
        %   not included in the path yet.
        %   warning ( newline ...
        %           + 'pm.sys.path.abs:NoMex', ...
        %           + ['pm.sys.path.abs: Using slow MATLAB-version instead of fast Mex.', ...
        %           + char(10), 'Compile: InstallMex pm.sys.path.abs.c'] ...
        %           + newline ...
        %           );
        %end
        % DATAREAD is deprecated in 2011b, but still available. In MATLAB 6.5, REGEXP
        % does not know the 'split' command, therefore DATAREAD is preferred.
        hasDataRead = ~isempty(which('dataread'));
    end

    if isempty(path) % Accept empty matrix as input.
        if ischar(path) || isnumeric(path)
            path = cd;
            return;
        else
            error   ( newline ...
                    + "A non-empty input argument `path` must be a string, character, or cell array of such values." + newline ...
                    + newline ...
                    );
        end
    end

    if ~ischar(path)
        % Non-empty inputs must be strings
        error   ( newline ...
                + "A non-empty input argument `path` must be a string, character, or cell array of such values." + newline ...
                + newline ...
                );
    end

    if isWindows % Windows: --------------------------------------------------------

        dirsep = '\';
        path = strrep(path, '/', dirsep);

        % Remove the magic key on demand, it is appended finally again.

        if strncmp(path, '\\?\', 4)
            if strncmpi(path, '\\?\UNC\', 8)
                % [BUGFIX] 19-May-2019, Thanks HHang Li, "path(7:..." -> "path(8:..."
                path = ['\', path(8 : length(path))];  % Two leading backslashes!
            else
                path = path(5 : length(path));
            end
        end

        isUNC = strncmp(path, '\\', 2);
        FileLen = length(path);
        if isUNC == 0 % path is not a UNC path

            % Leading file separator means relative to current drive or base folder.

            ThePath = cd;
            if path(1) == dirsep
                if strncmp(ThePath, '\\', 2)
                    % Current directory is a UNC path.
                    sepInd  = strfind(ThePath, '\');
                    ThePath = ThePath(1 : sepInd(4));
                else
                    % drive letter only.
                    ThePath = ThePath(1 : 3);
                end
            end

            if FileLen < 2 || path(2) ~= ':'
                % Does not start with drive letter.
                if ThePath(length(ThePath)) ~= dirsep
                    if path(1) ~= dirsep
                        path = [ThePath, dirsep, path];
                    else
                        % path starts with separator.
                        path = [ThePath, path];
                    end
                else
                    % Current path ends with separator.
                    if path(1) ~= dirsep
                        path = [ThePath, path];
                    else
                        % path starts with separator.
                        ThePath(length(ThePath)) = [];
                        path = [ThePath, path];
                    end
                end
            elseif FileLen == 2 && path(2) == ':'
                % "C:" current directory on C!
                % "C:" is the current directory on the C-disk, even if the current
                % directory is on another disk! This was ignored in MATLAB 6.5, but
                % modern versions considers this strange behaviour.
                if strncmpi(ThePath, path, 2)
                    path = ThePath;
                else
                    try
                        backCD = cd;
                        path   = cd(cd(path));
                        cd(backCD);
                    catch me
                        if exist(path, 'dir')
                            % No idea what could cause an error then!
                            rethrow(me);
                        else
                            % Reply "K:\" for not existing disk.
                            path = [path, dirsep];
                        end
                    end
                end
            end
        end

    else % Linux, MacOS: ---------------------------------------------------

        dirsep = '/';
        path = strrep(path, '\', dirsep);

        if strcmp(path, '~') || strncmp(path, '~/', 2)
            % Home directory.
            homedir = getenv('HOME');
            if ~isempty(homedir)
                path(1) = [];
                path    = [homedir, path];
            end
        elseif strncmpi(path, dirsep, 1) == 0
            % Append relative path to current folder.
            ThePath = cd;
            if ThePath(length(ThePath)) == dirsep
                path = [ThePath, path];
            else
                path = [ThePath, dirsep, path];
            end
        end

    end

    % Care for "\." and "\.." - no efficient algorithm, but the fast Mex is recommended at all!

    if ~isempty(strfind(path, [dirsep, '.']))

        if isWindows
            if strncmp(path, '\\', 2)
                % UNC path
                index = strfind(path, '\');
                if length(index) < 4
                    % UNC path without separator after the folder.
                    return;
                end
                drive = path(1:index(4));
                path(1:index(4)) = [];
            else
                drive = path(1:3);
                path(1:3) = [];
            end
        else
            % Unix, MacOS.
            isUNC = false;
            drive = dirsep;
            path(1) = [];
        end

        hasTrailFSep = (path(length(path)) == dirsep);
        if hasTrailFSep
            path(length(path)) = [];
        end

        if hasDataRead
            if isWindows
                % Need "\\" as separator.
                C = dataread('string', path, '%s', 'delimiter', '\\');  %#ok<REMFF1>
            else
                C = dataread('string', path, '%s', 'delimiter', dirsep);  %#ok<REMFF1>
            end
        else
            % Use the slower REGEXP, when DATAREAD is not available anymore.
            C = regexp(path, dirsep, 'split');
        end

        % Remove '\.\' directly without side effects.
        C(strcmp(C, '.')) = [];

        % Remove '\..' with the parent recursively.
        R = 1 : length(C);
        for dd = reshape(find(strcmp(C, '..')), 1, [])
            index    = find(R == dd);
            R(index) = [];
            if index > 1
                R(index - 1) = [];
            end
        end

        if isempty(R)
            path = drive;
            if isUNC && ~hasTrailFSep
                path(length(path)) = [];
            end
        elseif isWindows
            % If you have CStr2String, use the faster.
            % path = CStr2String(C(R), dirsep, hasTrailFSep);
            path = sprintf('%s\\', C{R});
            if hasTrailFSep
                path = [drive, path];
            else
                path = [drive, path(1:length(path) - 1)];
            end
        else
            % Unix.
            path = [drive, sprintf('%s/', C{R})];
            if ~hasTrailFSep
                path(length(path)) = [];
            end
        end

    end

    % "Very" long names under Windows.

    if isWindows
        if ~ischar(style)
            error   ( newline ...
                    + "The optional input argument `style` must be a string or character-valued." + newline ...
                    + newline ...
                    );
        end
        if (strncmpi(style, 'a', 1) && length(path) >= MAX_PATH) || ...
                strncmpi(style, 'f', 1)
            % Do not use [isUNC] here, because this concerns the input, which can
            % '.\path', while the current directory is an UNC path.
            if strncmp(path, '\\', 2)  % UNC path
                path = ['\\?\UNC', path(2:end)];
            else
                path = ['\\?\', path];
            end
        end
    end

end