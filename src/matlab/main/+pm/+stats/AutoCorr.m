%>  \brief
%>  This is the base class for generating objects containing
%>  information about autocorrelation of the input data.<br>
%>
%>  \note
%>  This is convenience class for easy computation
%>  of autocorrelation and its storage all in one place.<br>
%>  The primary advantage of this class over the MATLAB
%>  intrinsic functions is in the ability of this class
%>  to compute the result for input dataframe table and
%>  return the results always in MATLAB ``table`` format.<br>
%>  See the documentation of the class constructor below.<br>
%>
%>  \note
%>  See the documentation of the class constructor below.<br>
%>
%>  \final{AutoCorr}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:06 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef AutoCorr < pm.matlab.Handle

    properties(Access = public)
        %>
        %>  ``numlag``
        %>
        %>  The positive scalar MATLAB whole number representing the
        %>  number of lags for which the autocorrelation must be computed.<br>
        %>  This argument is directly passed to the corresponding argument
        %>  of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.<br>
        %>
        numlag = [];
        %>
        %>  ``numstd``
        %>
        %>  The positive scalar MATLAB double representing the
        %>  number of standard deviations to be used in computing the
        %>  lower and upper significance levels of the autocorrelation.<br>
        %>  This argument is directly passed to the corresponding argument
        %>  of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.<br>
        %>
        numstd = 1;
        %>
        %>  ``bnd``
        %>
        %>  The MATLAB table of rank ``2`` serving as a convenient
        %>  storage component each column of which corresponds to the
        %>  absolute ``numstd``-significance bound for the corresponding
        %>  computed autocorrelations of the input data (that is optionally
        %>  stored in the ``val`` component of the parent object).<br>
        %>  The values represent the ``numstd``-sigma standard
        %>  errors on the computed autocorrelations.<br>
        %>  Any autocorrelation value bound within these limits can be
        %>  considered random fluctuations at ``numstd``-sigma confidence level.<br>
        %>  This component is automatically populated when constructing
        %>  an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  It must be populated manually at all other times.<br>
        %>
        bnd = [];
        %>
        %>  ``lag``
        %>
        %>  The MATLAB vector integers representing the
        %>  lags for which the autocorrelation is computed.<br>
        %>  This component is automatically populated when constructing
        %>  an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  It must be populated manually at all other times.<br>
        %>
        lag = [];
        %>
        %>  ``val``
        %>
        %>  The MATLAB table of rank ``2`` serving as a
        %>  convenient storage component for the autocorrelation.<br>
        %>  This component is automatically populated when constructing
        %>  an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  It must be populated manually at all other times.<br>
        %>
        val = [];
    end

    methods(Access = public)

        %>  \brief
        %>  Return an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  This is the constructor of the [pm.stats.AutoCorr](@ref AutoCorr) class.<br>
        %>
        %>  \param[in]  df      :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose autocorrelation must be computed.<br>
        %>                          (**optional**. If missing, the autocorrelation will not be computed.)
        %>  \param[in]  numlag  :   The input scalar MATLAB positive whole number representing the number
        %>                          of lags for which the autocorrelation must be computed.<br>
        %>                          This argument is directly passed to the corresponding argument
        %>                          of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.<br>
        %>  \param[in]  numstd  :   The input positive scalar MATLAB double representing the
        %>                          number of standard deviations to be used in computing the
        %>                          lower and upper significance levels of the autocorrelation.<br>
        %>                          This argument is directly passed to the corresponding argument
        %>                          of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.<br>
        %>                          (**optional**, default = ``1``)
        %>
        %>  \return
        %>  ``self``            :   The output object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>
        %>  \interface{AutoCorr}
        %>  \code{.m}
        %>
        %>      acf = pm.stats.AutoCorr()
        %>      acf = pm.stats.AutoCorr(df)
        %>      acf = pm.stats.AutoCorr(df, numlag)
        %>      acf = pm.stats.AutoCorr(df, [], numstd)
        %>      acf = pm.stats.AutoCorr(df, numlag, numstd)
        %>
        %>  \endcode
        %>
        %>  \example{AutoCorr}
        %>  \include{lineno} example/stats/AutoCorr/main.m
        %>  \vis{AutoCorr}
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.png width=700
        %>
        %>  \final{AutoCorr}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:10 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = AutoCorr(df, numlag, numstd)
            if nargin < 3
                numstd = 1;
            end
            if nargin < 2
                numlag = [];
            end
            if 0 < nargin
                if  isempty(numlag)
                    numlag = size(df, 1) - 1;
                end
                [self.val, self.lag, self.bnd] = self.get(df, numlag, numstd);
            else
                self.numlag = numlag;
                self.numstd = numstd;
            end
        end

        %>  \brief
        %>  Return the autocorrelation of the input
        %>  data from a lag of ``0`` to ``numlag``.<br>
        %>
        %>  \details
        %>  This is a dynamic method of the [pm.stats.AutoCorr](@ref AutoCorr) class.<br>
        %>  This method automatically stores any input information
        %>  in the corresponding components of the parent object.<br>
        %>  However, any components of the parent object
        %>  corresponding to the output of this method
        %>  must be set explicitly manually.<br>
        %>
        %>  \param[inout]   self    :   The input/output parent object of class [pm.stats.AutoCorr](@ref AutoCorr)
        %>                              which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      df      :   The input MATLAB matrix or table of rank ``2``
        %>                              containing the data as ``ncol`` columns of ``nrow``
        %>                              observations whose autocorrelation must be computed.<br>
        %>  \param[in]      numlag  :   The input positive scalar MATLAB integer representing the
        %>                              number of lags to be used in computing the autocorrelation.<br>
        %>                              The default value will be used if the input ``numlag``
        %>                              is unspecified or empty ``[]``.<br>
        %>                              (**optional**, default = ``size(df, 1) - 1``)
        %>  \param[in]      numstd  :   The input positive scalar MATLAB double representing the
        %>                              number of standard deviations to be used in computing the
        %>                              lower and upper significance levels of the autocorrelation.<br>
        %>                              This argument is directly passed to the corresponding argument
        %>                              of the MATLAB intrinsic ``autocorr()`` in the Econometrics Toolbox.<br>
        %>                              The default value will be used if the input ``numstd`` is empty ``[]``.<br>
        %>                              (**optional**, default = ``1``)
        %>
        %>  \return
        %>  ``val``                 :   The output MATLAB ``table`` of size ``numlag + 1``
        %>                              containing the autocorrelation from lag ``0`` to ``numlag``.<br>
        %>  ``lag``                 :   The output MATLAB ``table`` of size ``numlag + 1``
        %>                              containing the autocorrelation lags from ``0`` to ``numlag``.<br>
        %>  ``bnd``                 :   The output MATLAB ``table`` of ``size(df, 2)`` rows by one column
        %>                              containing the absolute ``numstd``-significance level of the
        %>                              computed autocorrelations. Any autocorrelation value whose
        %>                              magnitude is smaller than the corresponding ``bnd`` element
        %>                              can be considered insignificant and mere fluctuation.<br>
        %>
        %>  \interface{get}
        %>  \code{.m}
        %>
        %>      acf = pm.stats.AutoCorr()
        %>      [acf.val, acf.lag, acf.bnd] = acf.get(df)
        %>      [acf.val, acf.lag, acf.bnd] = acf.get(df, [])
        %>      [acf.val, acf.lag, acf.bnd] = acf.get(df, numlag)
        %>      [acf.val, acf.lag, acf.bnd] = acf.get(df, [], numstd)
        %>      [acf.val, acf.lag, acf.bnd] = acf.get(df, numlag, numstd)
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the documentation of the class constructor
        %>  [pm.stats.AutoCorr](@ref AutoCorr::AutoCorr) for example usage.<br>
        %>
        %>  \final{get}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:14 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function [val, lag, bnd] = get(self, df, numlag, numstd)
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
                numlag = [];
            end
            self.numlag = numlag;
            self.numstd = numstd;
            if isa(df, "table")
                data = table2array(df);
            else
                data = df;
            end
            if ~isempty(self.numlag)
                numlag = self.numlag;
            else
                numlag = size(data, 1) - 1;
            end
            val = zeros(numlag + 1, size(data, 2));
            lag = zeros(numlag + 1, 1);
            bnd = zeros(size(data, 2), 1);
            try
                for icol = 1 : size(data, 2)
                    [val(:, icol), lag, bounds] = autocorr(data(:, icol), numlag);
                    bnd(icol) = abs(bounds(1)) * numstd;
                end
                val = array2table(val);
            catch me
                %val = NaN(numlag, size(data, 2));
                warning ( newline ...
                        + string(me.identifier) + " : " + string(me.message) + newline ...
                        + "skipping the autocorrelation computation..." + newline ...
                        + newline ...
                        );
                return;
            end
            if isa(df, "table")
                val.Properties.VariableNames = df.Properties.VariableNames;
            end
        end

    end

end