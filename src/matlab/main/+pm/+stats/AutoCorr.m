%
%   This is the base class for generating objects containing
%   information about autocorrelation of the input data.
%
%   \note
%
%       This is convenience class for easy computation
%       of autocorrelation and its storage all in one place.
%       The primary advantage of this class over the MATLAB
%       intrinsic functions is in the ability of this class
%       to compute the result for input dataframe table and
%       return the results always in MATLAB ``table`` format.
%
%       See the documentation of the class constructor below.
%
%   \note
%       See the documentation of the class constructor below.
%
%>  \return
%       See the documentation of the class constructor below.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       acf = pm.stats.AutoCorr()
%       acf = pm.stats.AutoCorr(df)
%       acf = pm.stats.AutoCorr(df, numlags)
%       acf = pm.stats.AutoCorr(df, [], numstd)
%       acf = pm.stats.AutoCorr(df, numlags, numstd)
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef AutoCorr < pm.matlab.Handle

    properties(Access = public)
        %
        %       numstd
        %
        %           The positive scalar MATLAB double representing the
        %           number of lags for which the autocorrelation must be computed.
        %           This argument is directly passed to the corresponding argument
        %           of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.
        %
        numlags = [];
        %
        %       numstd
        %
        %           The positive scalar MATLAB double representing the
        %           number of standard deviations to be used in computing the
        %           lower and upper significance levels of the autocorrelation.
        %           This argument is directly passed to the corresponding argument
        %           of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.
        %
        numstd = 1;
        %
        %   bnd
        %
        %       The MATLAB table of rank ``2`` serving as a convenient
        %       storage component each column of which corresponds to the
        %       absolute ``numstd``-significance bound for the corresponding
        %       computed autocorrelations of the input data (that is optionally
        %       stored in the ``val`` component of the parent object).
        %       The values represent the ``numstd``-sigma standard
        %       errors on the computed autocorrelations.
        %       Any autocorrelation value bound within these limits can be
        %       considered random fluctuations at ``numstd``-sigma confidence level.
        %       This component is automatically populated at the time
        %       of constructing an object of class ``pm.stats.AutoCorr``.
        %       It must be populated manually at all other times.
        %
        bnd = [];
        %
        %   lag
        %
        %       The MATLAB vector integers representing the
        %       lags for which the autocorrelation is computed.
        %       This component is automatically populated at the time
        %       of constructing an object of class ``pm.stats.AutoCorr``.
        %       It must be populated manually at all other times.
        %
        lag = [];
        %
        %   val
        %
        %       The MATLAB table of rank ``2`` serving as a
        %       convenient storage component for the autocorrelation.
        %       This component is automatically populated at the time
        %       of constructing an object of class ``pm.stats.AutoCorr``.
        %       It must be populated manually at all other times.
        %
        val = [];
    end

    methods(Access=public)

        %
        %   Return an object of class ``pm.stats.AutoCorr``.
        %
        %   This is the constructor of the ``pm.stats.AutoCorr`` class.
        %
        %       df
        %
        %           The input MATLAB matrix or table of rank ``2``
        %           containing the data as ``ncol`` columns of ``nrow``
        %           observations whose autocorrelation must be computed.
        %           (**optional**. If missing, the autocorrelation will not be computed.)
        %
        %       numstd
        %
        %           The input positive scalar MATLAB double representing the
        %           number of standard deviations to be used in computing the
        %           lower and upper significance levels of the autocorrelation.
        %           This argument is directly passed to the corresponding argument
        %           of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.
        %
        %           (**optional**, default = ``1``)
        %>
        %>  \return
        %>
        %       self
        %
        %           The output object of class ``pm.stats.AutoCorr``.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       acf = pm.stats.AutoCorr()
        %       acf = pm.stats.AutoCorr(df)
        %       acf = pm.stats.AutoCorr(df, numlags)
        %       acf = pm.stats.AutoCorr(df, [], numstd)
        %       acf = pm.stats.AutoCorr(df, numlags, numstd)
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = AutoCorr(df, numlags, numstd)
            if nargin < 3
                numstd = 1;
            end
            if nargin < 2
                numlags = [];
            end
            if 0 < nargin
                if  isempty(numlags)
                    numlags = size(df, 1) - 1;
                end
                [self.val, self.lag, self.bnd] = self.get(df, numlags, numstd);
            else
                self.numlags = numlags;
                self.numstd = numstd;
            end
        end

        %
        %   Return the autocorrelation of the input
        %   data from a lag of ``0`` to ``numlags``.
        %
        %   This is a dynamic method of the ``pm.stats.AutoCorr`` class.
        %   This method automatically stores any input information
        %   in the corresponding components of the parent object.
        %   However, any components of the parent object
        %   corresponding to the output of this method
        %   must be set explicitly manually.
        %
        %       df
        %
        %           The input MATLAB matrix or table of rank ``2``
        %           containing the data as ``ncol`` columns of ``nrow``
        %           observations whose autocorrelation must be computed.
        %
        %       numlags
        %
        %           The input positive scalar MATLAB integer representing the
        %           number of lags to be used in computing the autocorrelation.
        %           The default value will be used if the input ``numlags``
        %           is unspecified or empty ``[]``.
        %
        %           (**optional**, default = ``size(df, 1) - 1``)
        %
        %       numstd
        %
        %           The input positive scalar MATLAB double representing the
        %           number of standard deviations to be used in computing the
        %           lower and upper significance levels of the autocorrelation.
        %           This argument is directly passed to the corresponding argument
        %           of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.
        %           The default value will be used if the input ``numstd`` is empty ``[]``.
        %
        %           (**optional**, default = ``1``)
        %>
        %>  \return
        %>
        %       val
        %
        %           The output MATLAB ``table`` of size ``numlags + 1``
        %           containing the autocorrelation from lag ``0`` to ``numlags``.
        %
        %       lag
        %
        %           The output MATLAB ``table`` of size ``numlags + 1``
        %           containing the autocorrelation lags from ``0`` to ``numlags``.
        %
        %       bnd
        %
        %           The output MATLAB ``table`` of ``size(df, 2)`` rows by one column
        %           containing the absolute ``numstd``-significance level of the
        %           computed autocorrelations. Any autocorrelation value whose
        %           magnitude is smaller than the corresponding ``bnd`` element
        %           can be considered insignificant and mere fluctuation.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       acf = pm.stats.AutoCorr()
        %       [acf.val, acf.lag, acf.bnd] = acf.get(df)
        %       [acf.val, acf.lag, acf.bnd] = acf.get(df, [])
        %       [acf.val, acf.lag, acf.bnd] = acf.get(df, numlags)
        %       [acf.val, acf.lag, acf.bnd] = acf.get(df, [], numstd)
        %       [acf.val, acf.lag, acf.bnd] = acf.get(df, numlags, numstd)
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function [val, lag, bnd] = get(self, df, numlags, numstd)
            if nargin < 2
                help("pm.stats.AutoCorr");
                error   ( newline ...
                        + "The input ``df`` argument is required for computing the autocorrelation." + newline ...
                        + newline ...
                        );
            end
            if  nargin < 4
                numstd = 1;
            end
            if  nargin < 3
                numlags = [];
            end
            self.numlags = numlags;
            self.numstd = numstd;
            if isa(df, "table")
                data = table2array(df);
            else
                data = df;
            end
            if ~isempty(self.numlags)
                numlags = self.numlags;
            else
                numlags = size(data, 1) - 1;
            end
            val = zeros(numlags + 1, size(data, 2));
            lag = zeros(numlags + 1, 1);
            bnd = zeros(size(data, 2), 1);
            try
                for icol = 1 : size(data, 2)
                    [val(:, icol), lag, bounds] = autocorr(data(:, icol), numlags);
                    bnd(icol) = abs(bounds(1)) * numstd;
                end
                val = array2table(val);
            catch ME
                %val = NaN(numlags, size(data, 2));
                warning ( newline ...
                        + string(ME.message) + newline ...
                        + newline ...
                        + "skipping the autocorrelation computation..." + newline ...
                        );
                return;
            end
            if isa(df, "table")
                val.Properties.VariableNames = df.Properties.VariableNames;
            end
        end

    end

end