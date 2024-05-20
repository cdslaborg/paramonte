%
%   This is the base class for generating objects with methods
%   and storage components for computing and storing the
%   correlation matrix of an input data.
%
%   \note
%
%       This is convenience class for easy computation
%       of correlation and its storage all in one place.
%       The primary advantage of this class over the MATLAB
%       intrinsic functions is in the ability of this class
%       to compute the result for input dataframe table and
%       return the results always in MATLAB ``table`` format.
%
%       See the documentation of the class constructor below.
%
%   Attributes
%   ----------
%
%       See the documentation of the class constructor below.
%
%>  \return
%       See the documentation of the class constructor below.
%
%   Interface
%   ---------
%
%       mat = pm.stats.Cor()
%       mat = pm.stats.Cor(df)
%       mat = pm.stats.Cor(df, method)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Cor < pm.matlab.Handle

    properties(Access = public)
        %
        %   method
        %
        %       The scalar MATLAB string containing the
        %       method of computing the correlation matrix.
        %       It can be either:
        %
        %           "pearson"   : for computing the Pearson correlation matrix of the input data.
        %           "kendall"   : for computing the kendall rank correlation matrix of the input data.
        %           "spearman"  : for computing the Spearman rank correlation matrix of the input data.
        %
        %
        method = "pearson";
        %
        %   val
        %
        %       The MATLAB table of rank ``2`` serving as a
        %       convenient storage component for the correlation matrix.
        %       This component is automatically populated at the time of
        %       constructing an object of class ``pm.stats.Cor``.
        %       It must be populated manually at all other times.
        %
        val = [];
    end

    methods(Access=public)

        %
        %   Return an object of class ``Cor``.
        %
        %   This is the constructor of the ``Cor`` class.
        %
        %   Parameters
        %   ----------
        %
        %       df
        %
        %           The input MATLAB matrix or table of rank ``2``
        %           containing the data as ``ncol`` columns of ``nrow``
        %           observations whose correlation matrix must be computed.
        %           (**optional**. If missing, the correlation matrix will not be computed.)
        %
        %       method
        %
        %           The input scalar MATLAB string that can be either:
        %
        %               "pearson"   : for computing the Pearson correlation matrix of the input data.
        %               "kendall"   : for computing the kendall rank correlation matrix of the input data.
        %               "spearman"  : for computing the Spearman rank correlation matrix of the input data.
        %
        %           (**optional**, default = ``"pearson"``)
        %
        %   Returns
        %   -------
        %
        %       self
        %
        %           The output object of class ``pm.stats.Cor``.
        %
        %   Interface
        %   ---------
        %
        %       mat = pm.stats.Cor()
        %       mat = pm.stats.Cor(df)
        %       mat = pm.stats.Cor(df, method)
        %
        %   LICENSE
        %   -------
        %
        %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
        %
        function self = Cor(df, method)
            if nargin == 2
                self.method = method;
            end
            if 0 < nargin
                self.val = self.get(df, self.method);
            end
        end

        %
        %   Return the correlation matrix of the input data.
        %
        %   This is a dynamic method of the ``Cor`` class.
        %   This method automatically stores any input information
        %   in the corresponding components of the parent object.
        %   However, any components of the parent object
        %   corresponding to the output of this method
        %   must be set explicitly manually.
        %
        %   Parameters
        %   ----------
        %
        %       df
        %
        %           The input MATLAB matrix or table of rank ``2``
        %           containing the data as ``ncol`` columns of ``nrow``
        %           observations whose correlation matrix must be computed.
        %
        %       method
        %
        %           The input scalar MATLAB string that can be either:
        %
        %               "pearson"   : for computing the Pearson correlation matrix of the input data.
        %               "kendall"   : for computing the kendall rank correlation matrix of the input data.
        %               "spearman"  : for computing the Spearman rank correlation matrix of the input data.
        %
        %           (**optional**, default = ``pm.stats.Cor.method``)
        %
        %   Returns
        %   -------
        %
        %       val
        %
        %           The output MATLAB ``table`` containing the correlation matrix.
        %
        %   Interface
        %   ---------
        %
        %       mat = pm.stats.Cor()
        %       mat.val = mat.get(df)
        %       mat.val = mat.get(df, method)
        %
        %   LICENSE
        %   -------
        %
        %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
        %
        function val = get(self, df, method)
            if nargin < 2
                help("pm.stats.Cor");
                error   ( newline ...
                        + "The input ``df`` argument is required for computing the correlation matrix." + newline ...
                        + newline ...
                        );
            end
            if nargin < 3
                method = "pearson";
            end
            self.method = method;
            if isa(df, "table")
                data = table2array(df);
            else
                data = df;
            end
            try
                val = array2table(corr(data, "type", self.method));
            catch ME
                val = NaN(size(data, 2), size(data, 2));
                warning ( newline ...
                        + string(ME.message) + newline ...
                        + newline ...
                        + "skipping the correlation matrix computation..." + newline ...
                        );
                return;
            end
            if isa(df, "table")
                val.Properties.VariableNames = df.Properties.VariableNames;
            end
        end

    end

end