%>  \brief
%>  Expand wildcards for files and directory names.
%>
%>  \details
%>  Pattern matching of file and directory names, based on wildcard
%>  characters. This function is similar to wildcard expansion performed
%>  by the Unix shell and Python glob.glob function, but it can handle
%>  more types of wildcards.<br>
%>  The following list highlights the key differences
%>  between this function and the MATLAB intrinsic ``dir()``.<br>
%>      -   ``glob()`` supports wildcards for directories.
%>      -   ``glob()`` returns the directory part of ``pattern``.
%>      -   ``glob()`` returns a cell array of matching names.
%>      -   ``glob()`` does not return hidden files and directories that
%>          start with ``'.'`` unless explicitly specified in ``pattern``.
%>      -   ``glob()`` does not return ``'.'`` and ``'..'`` unless explicitly specified in ``pattern``.
%>      -   ``glob()`` adds a trailing file separator to directory names.
%>      -   ``glob()`` does not return the contents of a directory when a directory
%>          is specified. To return contents of a directory, add a trailing ``'/*'``.
%>      -   ``glob()`` returns only directory names when a trailing file separator is specified.
%>      -   On Windows, ``glob()`` is not case sensitive, but it returns matching
%>          names exactly in the case as they are defined on the filesystem.
%>          Case of host and sharename of a UNC path and case of drive
%>          letters will be returned as specified in ``pattern``.
%>
%>  \param[in]  pattern :   The input scalar MATLAB string, containing the search pattern.
%>                          Wildcards may be used for basenames and for the directory parts.
%>                          If pattern contains directory parts, then these will be included in the output ``pathList``.<br>
%>                          Following wildcards can be used:<br>
%>                              *        match zero or more characters
%>                              ?        match any single character
%>                              [ab12]   match one of the specified characters
%>                              [^ab12]  match none of the specified characters
%>                              [a-z]    match one character in range of characters
%>                              {a,b,c}  matches any one of strings a, b or c<br>
%>                              All above wildcards do not match a file separator.<br>
%>                              **       match zero or more characters including file separators.
%>                                       This can be used to match zero or more directory parts
%>                                       and will recursively list matching names.
%>  
%>                              \warning
%>                              Symbolic linked directories or junctions may
%>                              cause an infinite loop when using the ``**``.
%>  
%>  \param[in]  anycase :   The input scalar MATLAB logical.
%>                          If ``true``, the search will be case-sensitive.
%>                          If ``false``, the search will be case-insensitive.
%>                          On Windows, ``anycase`` is always reset to ``true`` even if user-specified.
%>                          (**optional**. default = ``false`` on Unix and ``true`` on Windows.)
%>  
%>  \return
%>  `pathList`          :   The output MATLAB cell array of strings containing the files
%>                          or directories that match the path specified by string ``pattern``.<br>
%>  `isdirList`         :   The output MATLAB cell array of the same size as ``pathList``,
%>                          each element of which is a MATLAB logical value that is ``true`` if
%>                          and only if the corresponding element of ``pathList`` is a directory.
%>
%>  \interface{glob}
%>  \code{.m}
%>
%>      [pathList, isdirList] = pm.sys.path.glob(pattern)
%>      [pathList, isdirList] = pm.sys.path.glob(pattern, anycase)
%>
%>  \endcode
%>
%>  \example{getBorder}
%>
%>      glob *.m        list all .m files in current directory.
%>
%>      glob baz/*      list all files and directories in subdirectory 'baz'.
%>
%>      glob b*/*.m     list all .m files in subdirectory names starting
%>                      with 'b'. The list will include the names of the
%>                      matching subdirectories.
%>
%>      glob ?z*.m      list all .m files where the second character
%>                      is 'z'.
%>
%>      glob baz.[ch]   matches baz.c and baz.h
%>
%>      glob test.[^ch] matches test.a but not test.c or test.h
%>
%>      glob demo.[a-c] matches demo.a, demo.b, and demo.c
%>
%>      glob test.{foo,bar,baz} matches test.foo, test.bar, and test.baz
%>
%>      glob .*         list all hidden files in current directory,
%>                      excluding '.' and '..'
%>
%>      glob */         list all subdirectories.
%>
%>      glob **         recursively list all files and directories,
%>                      starting in current directory (current directory
%>                      name, hidden files and hidden directories are
%>                      excluded).
%>
%>      glob **.m       list all m-files anywhere in directory tree,
%>                      including m-files in current directory. This
%>                      is equivalent with '**/*.m'.
%>
%>      glob foo/**/    recursively list all directories, starting in
%>                      directory 'foo'.
%>
%>      glob **/.svn/   list all .svn directories in directory tree.
%>
%>      glob **/.*/**   recursively list all files in hidden directories
%>                      only.
%>
%>      [r,d]=glob('**')
%>      r(~d)           get all files in directory tree.
%>
%>  \final{glob}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:24 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>      Copyright (c) 2013, Peter van den Biggelaar
%>      All rights reserved.
%>
%>      Redistribution and use in source and binary forms, with or without
%>      modification, are permitted provided that the following conditions are
%>      met:
%>
%>          * Redistributions of source code must retain the above copyright
%>            notice, this list of conditions and the following disclaimer.
%>          * Redistributions in binary form must reproduce the above copyright
%>            notice, this list of conditions and the following disclaimer in
%>            the documentation and/or other materials provided with the distribution
%>
%>      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%>      AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%>      IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%>      ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%>      LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%>      CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%>      SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%>      INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%>      CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%>      ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%>      POSSIBILITY OF SUCH DAMAGE.
%>
function [pathList, isdirList] = glob(pattern, anycase)
    if isstring(pattern)
        pattern = convertStringsToChars(pattern);
    end
    % check pattern input
    if ischar(pattern)
        if isempty(pattern)
            % return when pattern is empty
            pathList = cell(0);
            isdirList = false(0);
            return
        elseif size(pattern,1)>1
            error('glob:invalidInput', 'pattern must be a single string.')
        end
    else
        error('glob:invalidInput', 'pattern must be a string.')
    end
    % check anycase option
    if nargin == 2
        pm.introspection.verify(anycase, "logical", 1, "anycase");
    else
        % Windows is not case sensitive
        % Unix is case sensitive
        anycase = ispc;
    end
    % define function handle to regular expression function for the specified case sensitivity
    if anycase
        regexp_fhandle = @regexpi;
    else
        regexp_fhandle = @regexp;
    end
    % only use forward slashes as file separator to prevent escaping backslashes in regular expressions
    filespec = strrep(pattern, '\', '/');
    % split pathroot part from pattern
    if strncmp(filespec, '//',2)
        if ispc
            % pattern specifies a UNC path
            % It is not allowed to get a directory listing of share names of a
            % host with the DIR command.
            % pathroot will contains e.g. //host/share/
            pathroot = regexprep(filespec, '(^//+[^/]+/[^/]+/)(.*)', '$1');
            filespec = regexprep(filespec, '(^//+[^/]+/[^/]+/)(.*)', '$2');
        else
            % for Unix, multiple leading file separators are equivalent with a single file separator
            filespec = regexprep(filespec, '^/*', '/');
        end
    elseif strncmp(filespec, '/', 1)
        % pattern specifies a absolute path
        pathroot = '/';
        filespec(1) = [];
    elseif ispc && numel(filespec)>=2 && filespec(2) == ':'
        % pattern specifies a absolute path starting with a drive letter
        % check for a fileseparator after ':'. e.g. 'C:\'
        if numel(filespec)<3 || filespec(3)~='/'
            error('glob:invalidInput','Drive letter must be followed by '':\''.')
        end
        pathroot = filespec(1:3);
        filespec(1:3) = [];
    else
        % pattern specifies a relative path
        pathroot = './';
    end
    % replace multiple file separators by a single file separator
    filespec = regexprep(filespec, '/+', '/');
    % replace 'a**' with 'a*/**', where 'a' can be any character but not '/'
    filespec = regexprep(filespec, '([^/])(\.\*\.\*)', '$1\*/$2');
    % replace '**a' with '**/*a', where a can be any character but not '/'
    filespec = regexprep(filespec, '(\.\*\.\*)([^/])', '$1/\*$2');
    % split filespec into chunks at file separator
    chunks = strread(filespec, '%s', 'delimiter', '/'); %#ok<FPARK>
    % add empty chunk at the end when filespec ends with a file separator
    if ~isempty(filespec) && filespec(end)=='/'
        chunks{end+1} = '';
    end
    % translate chunks to regular expressions
    for i=1:numel(chunks)
        chunks{i} = glob2regexp(chunks{i});
    end
    % determine file list using LS_REGEXP
    % this function requires that PATHROOT does not to contain any wildcards
    if ~isempty(chunks)
        list = ls_regexp(regexp_fhandle, pathroot, chunks{1:end});
    else
        list = {pathroot};
    end
    if strcmp(pathroot, './')
        % remove relative pathroot from result
        list = regexprep(list, '^\./', '');
    end
    if nargout == 2
        % determine directories by checking for '/' at the end
        I = regexp(list', '/$');
        isdirList = ~cellfun('isempty', I);
    end
    % convert to standard file separators for PC
    if ispc
        list = strrep(list, '/', '\');
    end
    % return output
    if nargout == 0
        if ~isempty(list)
            % display list
            disp(string(list))
        else
            disp(['''' pattern ''' not found.']);
        end
    else
        pathList = string(list');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function regexp_str = glob2regexp(glob_str)
        % translate glob_str to regular expression string
        % initialize
        regexp_str  = '';
        in_curlies  = 0;        % is > 0 within curly braces
        % handle characters in glob_str one-by-one
        for c = glob_str

            if any(c=='.()|+^$@%')
                % escape simple special characters
                regexp_str = [regexp_str '\' c]; %#ok<AGROW>

            elseif c=='*'
                % '*' should not match '/'
                regexp_str = [regexp_str '[^/]*']; %#ok<AGROW>

            elseif c=='?'
                % '?' should not match '/'
                regexp_str = [regexp_str '[^/]']; %#ok<AGROW>

            elseif c=='{'
                regexp_str = [regexp_str '(']; %#ok<AGROW>
                in_curlies = in_curlies+1;
            elseif c=='}' && in_curlies
                regexp_str = [regexp_str ')']; %#ok<AGROW>
                in_curlies = in_curlies-1;
            elseif c==',' && in_curlies
                regexp_str = [regexp_str '|']; %#ok<AGROW>

            else
                regexp_str = [regexp_str c]; %#ok<AGROW>
            end
        end
        % replace original '**' (that has now become '[^/]*[^/]*') with '.*.*'
        regexp_str = strrep(regexp_str, '[^/]*[^/]*', '.*.*');
    end

    function L = ls_regexp(regexp_fhandle, path, varargin)
        % List files that match PATH/r1/r2/r3/... where PATH is a string without
        % any wildcards and r1..rn are regular expresions that contain the parts of
        % a filespec between the file separators.
        % L is a cell array with matching file or directory names.
        % REGEXP_FHANDLE contain a file handle to REGEXP or REGEXPI depending
        % on specified case sensitivity.

        % if first regular expressions contains '**', examine complete file tree
        if nargin >= 3 && any(regexp(varargin{1}, '\.\*\.\*'))
            L = ls_regexp_tree(regexp_fhandle, path, varargin{:});

        else
            % get contents of path
            list = dir(path);

            if nargin >= 3
                if strcmp(varargin{1}, '\.') || strcmp(varargin{1}, '\.\.')
                    % keep explicitly specified '.' or '..' in first regular expression
                    if ispc && ~any(strcmp({list.name}, '.'))
                        % fix strange windows behaviour: root of a volume has no '.' and '..'
                        list(end+1).name = '.';
                        list(end).isdir = true;
                        list(end+1).name = '..';
                        list(end).isdir = true;
                    end
                else
                    % remove '.' and '..'
                    list(strcmp({list.name},'.')) = [];
                    list(strcmp({list.name},'..')) = [];

                    % remove files starting with '.' specified in first regular expression
                    if ~strncmp(varargin{1},'\.',2)
                        % remove files starting with '.' from list
                        list(strncmp({list.name},'.',1))  = [];
                    end
                end
            end

            % define shortcuts
            list_isdir = [list.isdir];
            list_name = {list.name};

            L = {};  % initialize
            if nargin==2    % no regular expressions
                % return filename
                if ~isempty(list_name)
                    % add a trailing slash to directories
                    trailing_fsep = repmat({''}, size(list_name));
                    trailing_fsep(list_isdir) = {'/'};
                    L = strcat(path, list_name, trailing_fsep);
                end
            elseif nargin==3    % last regular expression
                % return list_name matching regular expression
                I = regexp_fhandle(list_name, ['^' varargin{1} '$']);
                I = ~cellfun('isempty', I);
                list_name = list_name(I);
                list_isdir = list_isdir(I);
                if ~isempty(list_name)
                    % add a trailing slash to directories
                    trailing_fsep = repmat({''}, size(list_name));
                    trailing_fsep(list_isdir) = {'/'};
                    L = strcat(path, list_name, trailing_fsep);
                end

            elseif nargin==4 && isempty(varargin{2})
                % only return directories when last regexp is empty
                % return list_name matching regular expression and that are directories
                I = regexp_fhandle(list_name, ['^' varargin{1} '$']);
                I = ~cellfun('isempty', I);
                % only return directories
                list_name = list_name(I);
                list_isdir = list_isdir(I);
                if any(list_isdir)
                    % add a trailing file separator
                    L = strcat(path, list_name(list_isdir), '/');
                end
            else
                % traverse for list_name matching regular expression
                I = regexp_fhandle(list_name, ['^' varargin{1} '$']);
                I = ~cellfun('isempty', I);
                for name = list_name(I)
                    L = [L   ls_regexp(regexp_fhandle, [path char(name) '/'], varargin{2:end})]; %#ok<AGROW>
                end
            end
        end
    end

    function L = ls_regexp_tree(regexp_fhandle, path, varargin)
        % use this function when first argument of varargin contains '**'
        % build list of complete directory tree
        % if any regexp starts with '\.', keep hidden files and directories
        I = regexp(varargin, '^\\\.');
        I = ~cellfun('isempty', I);
        keep_hidden = any(I);
        list = dir_recur(path, keep_hidden);
        L = {list.name};
        % make one regular expression of all individual regexps
        expression = [regexptranslate('escape',path) sprintf('%s/', varargin{1:end-1}) varargin{end}];
        % note that /**/ must also match zero directories
        % replace '/**/' with (/**/|/)
        expression = regexprep(expression, '/\.\*\.\*/', '(/\.\*\.\*/|/)');
        % return matching names
        if ~isempty(varargin{end})
            % determing matching names ignoring trailing '/'
            L_no_trailing_fsep = regexprep(L, '/$', '');
            I = regexp_fhandle(L_no_trailing_fsep, ['^' expression '$']);
        else
            % determing matching names including trailing '/'
            I = regexp_fhandle(L, ['^' expression '$']);
        end
        I = cellfun('isempty', I);
        L(I) = [];
    end

    function d = dir_recur(startdir, keep_hidden)
        % determine recursive directory contents
        % get directory contents
        d = dir(startdir);
        % remove hidden files
        if keep_hidden
            % only remove '.' and '..'
            d(strcmp({d.name},'.'))  = [];
            d(strcmp({d.name},'..')) = [];
        else
            % remove all hidden files and directories
            d(strncmp({d.name},'.',1)) = [];
        end
        if ~isempty(d)
            % add trailing fileseparator to directories
            trailing_fsep = repmat({''}, size(d));
            trailing_fsep([d.isdir]) = {'/'};
            % prefix startdir to name and postfix fileseparator for directories
            dname = strcat(startdir, {d.name}, trailing_fsep');
            [d(:).name] = deal(dname{:});
            % recurse into subdirectories
            for subd = {d([d.isdir]).name}
                d = [d; dir_recur(char(subd), keep_hidden)]; %#ok<AGROW>
            end
        end
    end
end