%>  \brief
%>  This is the base class for generating objects
%>  that contain the **tabular contents** of a given file.
%>
%>  \details
%>  This class is meant to be primarily internally used
%>  by the ParaMonte library routines (e.g., samplers).
%>
%>  Tabular file comprises a header line followed by a table of data.<br>
%>  If there are more than one header line, only the last such line is
%>  considered as header and the rest are discarded as text.<br>
%>
%>  \note
%>  See the documentation of the class constructor.<br>
%>
%>  \note
%>  See below for information on the attributes (properties).
%>
%>  \note
%>  See below for information on the methods.
%>
%>  \final{FileContentsTabular}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:08 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsTabular < pm.io.FileContents

    properties(Access = public)
        %>
        %>  sep     :   The scalar MATLAB string representing
        %>              the field separator used in the file.<br>
        %>              This property is either specified by
        %>              the user or is inferred from the
        %>              file to be read.
        %>
        sep = [];
        %>
        %>  ncol    :   The scalar MATLAB integer representing
        %>              the number of columns identified in the file.
        %>
        ncol = 0;
        %>
        %>  nrow    :   The scalar MATLAB integer representing
        %>              the number of rows identified in the file.
        %>
        nrow = 0;
        %>
        %>  df      :   The scalar MATLAB table containing the
        %>              contents of the file as a dataframe.
        %>
        df = [];
    end

    methods(Access = public)

        %>  \brief
        %>  Construct and return a scalar object of class [pm.io.FileContentsTabular](@ref FileContentsTabular).
        %>
        %>  \details
        %>  This is the constructor of the class [pm.io.FileContentsTabular](@ref FileContentsTabular).
        %>
        %>  \param[in]  file    :   See the documentation of the corresponding argument of [pm.io.FileContents](@ref FileContents).
        %>  \param[in]  silent  :   See the documentation of the corresponding argument of [pm.io.FileContents](@ref FileContents).
        %>  \param[in]  sep     :   The input scalar MATLAB string containing the field separator used in the file.<br>
        %>                          Use ``"\t"`` to represent the horizontal tab character.<br>
        %>                          (**optional**. The default is inferred from the contents of the specified input ``file``.)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.io.FileContentsTabular](@ref FileContentsTabular).
        %>
        %>  \interface{FileContentsTabular}
        %>  \code{.m}
        %>
        %>      contents = pm.io.FileContentsTabular(file)
        %>      contents = pm.io.FileContentsTabular(file, [])
        %>      contents = pm.io.FileContentsTabular(file, silent)
        %>      contents = pm.io.FileContentsTabular(file, [], [])
        %>      contents = pm.io.FileContentsTabular(file, [], sep)
        %>      contents = pm.io.FileContentsTabular(file, silent, [])
        %>      contents = pm.io.FileContentsTabular(file, silent, sep)
        %>
        %>  \endcode
        %>
        %>  \final{FileContentsTabular}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:11 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
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
            if  pm.array.len(sep) == 0
                sep = "";
            else
                sep = string(sep);
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
                self.sep = string(self.sep);
            elseif length(convertStringsToChars(sep)) == 1 || strcmpi(sep, "\t")
                [d, self.sep, ~] = importdata(self.file, convertStringsToChars(sep));
                self.sep = string(self.sep);
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
                        + "The default or the inferred choice of field separator (""" + sep + """) may be incorrect." + newline ...
                        + "Alternatively, the specified file contents may be corrupted and non-tabular." + newline ...
                        + "If the specified tile path is a weblink, ensure the weblink corrently points" + newline ...
                        + "to the target raw tabular ASCII file and not an HTML page." + newline ...
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