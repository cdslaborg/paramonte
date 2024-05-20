%
%   This is the base class for generating objects
%   that contain the tabular contents of a given file.
%
%   This class is meant to be primarily internally used
%   by the ParaMonte library routines (e.g., samplers).
%
%   Tabular file comprises a header
%   line followed by a table of data.
%   If there are more than one header line,
%   only the last such line is considered as
%   header and the rest are discarded as text.
%
%       See the documentation of the class constructor.
%
%   Attributes
%   ----------
%
%       See below for information on the attributes (properties).
%
%   Methods
%   -------
%
%       See below for information on the methods.
%
%>  \return
%       An object of class pm.io.FileContentsTabular.
%
%   Interface
%   ---------
%
%       contents = pm.io.FileContentsTabular(file)
%       contents = pm.io.FileContentsTabular(file, [])
%       contents = pm.io.FileContentsTabular(file, silent)
%       contents = pm.io.FileContentsTabular(file, [], [])
%       contents = pm.io.FileContentsTabular(file, [], sep)
%       contents = pm.io.FileContentsTabular(file, silent, [])
%       contents = pm.io.FileContentsTabular(file, silent, sep)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsTabular < pm.io.FileContents
    properties(Access = public)
        %
        %   sep
        %
        %       The scalar MATLAB string representing
        %       the field separator used in the file.
        %       This property is either specified by
        %       the user or is inferred from the
        %       file to be read.
        %
        sep = [];
        %
        %   ncol
        %
        %       The scalar MATLAB integer representing
        %       the number of columns identified in the file.
        %
        ncol = 0;
        %
        %   nrow
        %
        %       The scalar MATLAB integer representing
        %       the number of rows identified in the file.
        %
        nrow = 0;
        %
        %   df
        %
        %       The scalar MATLAB table containing the
        %       contents of the file as a dataframe.
        %
        df = [];
    end

    methods(Access = public)

        %
        %   Return a scalar object of class ``pm.io.FileContentsTabular``.
        %
        %   This is the constructor of the class ``pm.io.FileContentsTabular``.
        %
        %   Parameters
        %   ----------
        %
        %       file
        %
        %           See the documentation of the corresponding argument of ``pm.io.FileContents``.
        %
        %       silent
        %
        %           See the documentation of the corresponding argument of ``pm.io.FileContents``.
        %
        %       sep
        %
        %           The input scalar MATLAB string
        %           containing the field separator used in the file.
        %           (**optional**. The default is inferred from the contents of the specified input ``file``.)
        %
        %   Returns
        %   -------
        %
        %       self
        %
        %           The output scalar object of class ``pm.io.FileContentsTabular``.
        %
        %   Interface
        %   ---------
        %
        %       contents = pm.io.FileContentsTabular(file)
        %       contents = pm.io.FileContentsTabular(file, [])
        %       contents = pm.io.FileContentsTabular(file, silent)
        %       contents = pm.io.FileContentsTabular(file, [], [])
        %       contents = pm.io.FileContentsTabular(file, [], sep)
        %       contents = pm.io.FileContentsTabular(file, silent, [])
        %       contents = pm.io.FileContentsTabular(file, silent, sep)
        %
        %   LICENSE
        %   -------
        %
        %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
        %
        function self = FileContentsTabular(file, silent, sep)
            if  nargin < 3
                sep = [];
            end
            if  nargin < 2
                silent = [];
            end
            if  nargin < 1
                file = [];
            end
            self = self@pm.io.FileContents(file, silent);
            if ~isfile(self.file)
                error   ( newline ...
                        + "The specified input file:" + newline ...
                        + newline ...
                        + pm.io.tab + """" + self.file + """" + newline ...
                        + newline ...
                        + "does not exist." + newline ...
                        + newline ...
                        );
            end
            if  pm.array.len(sep) == 0
                [d, self.sep, ~] = importdata(self.file);
            elseif length(convertStringsToChars(sep)) == 1
                [d, self.sep, ~] = importdata(self.file, sep);
            else
                weblinks = pm.lib.weblinks();
                error   ( newline ...
                        + "Multicharacter separators are unsupported as of ParaMonte 2." + newline ...
                        + "You have specified sep = """ + sep + """." + newline ...
                        + "If this feature is desired, please report this request at:" + newline ...
                        + newline ...
                        + pm.io.tab + pm.web.href(weblinks.github.issues.url) + newline ...
                        + newline ...
                        );
            end
            if  isempty(d)
                error   ( newline ...
                        + "The specified input file is empty:" + newline ...
                        + newline ...
                        + pm.io.tab + self.file + newline ...
                        + newline ...
                        );
            elseif ~isfield(d, "colheaders") || ~isfield(d, "data")
                error   ( newline ...
                        + "The specified input file is not a table:" + newline ...
                        + newline ...
                        + pm.io.tab + self.file + newline ...
                        + newline ...
                        + "The choice of field separator (""" + sep + """) may be incorrect." + newline ...
                        + newline ...
                        );
            end
            self.ncol = size(d.data, 2);
            self.nrow = size(d.data, 1);
            self.df = array2table(d.data, 'VariableNames', d.colheaders);
            self.checkpoint();
            if ~self.silent
                disp("ncol = " + string(self.ncol) + ", nrow = " + string(self.nrow));
            end
        end

    end

end