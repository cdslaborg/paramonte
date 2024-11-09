%>  \brief
%>  This is the base class for generating objects containing
%>  information about autocorrelation of the input data and
%>  tools for visualizing it.<br>
%>
%>  \note
%>  This is convenience class for easy computation
%>  of autocorrelation and its storage all in one place.<br>
%>  The primary advantage of this class over the MATLAB
%>  intrinsic functions is in the ability of this class
%>  to compute the result for input dataframe table and
%>  return the results always in MATLAB ``table`` format.<br>
%>
%>  \note
%>  See the documentation of the class constructor
%>  [pm.stats.AutoCorr::AutoCorr](@ref AutoCorr::AutoCorr) below.<br>
%>
%>  \final{AutoCorr}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:06 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef AutoCorr < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``df``
        %>
        %>  A scalar object of class [pm.container.DataFrame](@ref DataFrame)
        %>  containing the user-specified data whose covariance must be computed.<br>
        %>
        df = [];
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
        %>  The MATLAB ``table`` of rank ``2`` serving as a convenient
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
        %>  This component is automatically populated when constructing
        %>  an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  It must be populated manually at all other times.<br>
        %>
        lag = [];
        %>
        %>  ``val``
        %>
        %>  The MATLAB ``table`` of rank ``2`` serving as a
        %>  convenient storage component for the autocorrelation.<br>
        %>  This component is automatically populated when constructing
        %>  an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  It must be populated manually at all other times.<br>
        %>
        %>  \note
        %>  The first column of ``vall`` always contains the set
        %>  of lags for which the autocorrelation is computed.<br>
        %>
        val = [];
        %>
        %>  ``vis``
        %>
        %>  The scalar MATLAB ``struct`` containing the set of
        %>  predefined visualizations for the output data.<br>
        %>
        vis = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Return an object of class [pm.stats.AutoCorr](@ref AutoCorr).<br>
        %>  This is the constructor of the [pm.stats.AutoCorr](@ref AutoCorr) class.<br>
        %>
        %>  \param[in]  dfref   :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose autocorrelation must be computed.<br>
        %>                          Ideally, the user would want to pass a reference
        %>                          to a dataframe (e.g., as a function handle ``@()df``)
        %>                          so that the data remains dynamically up-to-date.<br>
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
        %>      acf = pm.stats.AutoCorr([])
        %>      acf = pm.stats.AutoCorr(dfref)
        %>      acf = pm.stats.AutoCorr(dfref, numlag)
        %>      acf = pm.stats.AutoCorr(dfref, [], numstd)
        %>      acf = pm.stats.AutoCorr(dfref, numlag, numstd)
        %>
        %>  \endcode
        %>
        %>  \example{AutoCorr}
        %>  \include{lineno} example/stats/AutoCorr/main.m
        %>  \vis{AutoCorr}
        %>  <br><br>
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.plot.line.png width=700
        %>  <br><br>
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.tile.line.png width=700
        %>  <br><br>
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.cascade.line.1.png width=700
        %>  <br><br>
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.cascade.line.2.png width=700
        %>  <br><br>
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.cascade.line.3.png width=700
        %>  <br><br>
        %>  \image html example/stats/AutoCorr/AutoCorr.unifrnd.cascade.line.4.png width=700
        %>
        %>  \final{AutoCorr}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:10 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = AutoCorr(dfref, numlag, numstd)

            if  nargin < 3
                numstd = 1;
            end

            if  nargin < 2
                numlag = [];
            end

            if  nargin < 1
                self.numlag = numlag;
                self.numstd = numstd;
            else
                [self.val, self.bnd] = self.get(dfref, numlag, numstd);
            end

            self.setvis();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %>  \param[in]      dfref   :   The input MATLAB matrix or table of rank ``2``
        %>                              containing the data as ``ncol`` columns of ``nrow``
        %>                              observations whose autocorrelation must be computed.<br>
        %>                              Ideally, the user would want to pass a reference
        %>                              to a dataframe (e.g., as a function handle ``@()df``)
        %>                              so that the data remains dynamically up-to-date.<br>
        %>                              (**optional**. If missing, the contents of the corresponding
        %>                              internal component of the parent object will be used.)
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
        %>                              The first column of ``val`` always contains the set of ``numlag + 1``
        %>                              autocorrelation lags from ``0`` to ``numlag``.<br>
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
        %>      [acf.val, acf.bnd] = acf.get(dfref)
        %>      [acf.val, acf.bnd] = acf.get(dfref, [])
        %>      [acf.val, acf.bnd] = acf.get(dfref, numlag)
        %>      [acf.val, acf.bnd] = acf.get(dfref, [], numstd)
        %>      [acf.val, acf.bnd] = acf.get(dfref, numlag, numstd)
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
        function [val, bnd] = get(self, dfref, numlag, numstd)

            %%%%
            %%%% Define the data.
            %%%%

            if  nargin < 2
                dfref = [];
            end

            if ~isempty(dfref)
                self.df = pm.container.DataFrame(dfref);
            end

            dfcopy = self.df.copy();

            if ~isempty(dfcopy)
                data = dfcopy{:,:};
            else
                help("pm.stats.AutoCorr");
                error   ( newline ...
                        + "A non-mepty ``dfref`` attribute or input argument is required for computing the autocorrelation." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Define the number of lags.
            %%%%

            if  nargin < 3
                numlag = [];
            end

            if ~isempty(numlag)
                self.numlag = numlag;
            elseif isempty(self.numlag)
                self.numlag = size(data, 1) - 1;
            end

            %%%%
            %%%% Define the number of standard deviations for the bounds.
            %%%%

            if  nargin < 4
                numstd = [];
            end

            if ~isempty(numstd)
                self.numstd = numstd;
            elseif isempty(self.numstd)
                self.numstd = 1;
            end

            %%%%
            %%%% Compute the ACF.
            %%%%

            bnd = zeros(size(data, 2), 1);
            val = zeros(self.numlag + 1, size(data, 2) + 1);

            try

                for icol = 1 : size(data, 2)
                    [val(:, icol + 1), val(:, 1), bounds] = autocorr(data(:, icol), "numlags", self.numlag);
                    bnd(icol) = abs(bounds(1)) * self.numstd;
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

            val.Properties.VariableNames = ["Lag", "ACF(" + string(dfcopy.Properties.VariableNames) + ")"];

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Set up the visualization tools of the autocorrelation.<br>
        %>
        %>  \details
        %>  This is a dynamic ``Hidden`` method of the [pm.stats.AutoCorr](@ref AutoCorr) class.<br>
        %>  This method is inaccessible to the end users of the ParaMonte MATLAB library.<br>
        %>
        %>  \param[inout]   self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]      val     :   The input (reference of function handle returning a) MATLAB matrix or table of rank ``2``
        %>                              containing the computed autocorrelation to be visualized.<br>
        %>                              Ideally, the user would want to pass a reference to a dataframe (e.g., as a function handle ``@()df``)
        %>                              so that the data remains dynamically up-to-date.<br>
        %>                              (**optional**. If missing, the contents of the corresponding ``var`` attribute of the parent object will be used.)
        %>
        %>  \interface{get}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.AutoCorr(dfref, method)
        %>      mat.setvis(); % This method is automatically called within the object constructor.
        %>
        %>  \endcode
        %>
        %>  \final{get}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function setvis(self, val)

            if  1 < nargin
                self.val = val;
            end

            self.vis = struct();

            self.vis.plot = struct();
            self.vis.plot.line = pm.vis.PlotLine( @()self.val ...
                                                , "colx", 1, "coly", 2 : size(self.val, 2) ...
                                                , "legend", {"enabled", true} ..., "labels", self.val.Properties.VariableNames(2:end)} ...
                                                , "ylabel", {"txt", "Autocorrelation Function (ACF)"} ...
                                                , "xlabel", {"txt", "Autocorrelation Lag"} ...
                                                , "colormap", {"enabled", false} ...
                                                , "axes", {"xscale", "log"} ...
                                                , "ylim", [-1, 1] ...
                                                );

            self.vis.tile = struct();
            self.vis.tile.line = pm.vis.TileLine( @()self.val ...
                                                , "colx", 1, "coly", 2 : size(self.val, 2) ...
                                                , "ylabel", {"txt", "Autocorrelation Function (ACF)"} ...
                                                , "xlabel", {"txt", "Autocorrelation Lag"} ...
                                                , "colormap", {"enabled", false} ...
                                                , "legend", {"enabled", true} ...
                                                , "axes", {"xscale", "log"} ...
                                                , "tiledlayout", {"TileSpacing", "tight"} ...
                                                , "ylim", [-1, 1] ...
                                                );

            self.vis.cascade = struct();
            self.vis.cascade.line = pm.vis.CascadeLine  ( @()self.val ...
                                                        , "colx", 1, "coly", 2 : size(self.val, 2) ...
                                                        , "ylabel", {"txt", "Autocorrelation Function (ACF)"} ...
                                                        , "xlabel", {"txt", "Autocorrelation Lag"} ...
                                                        , "colormap", {"enabled", false} ...
                                                        , "legend", {"enabled", true} ...
                                                        , "axes", {"xscale", "log"} ...
                                                        , "ylim", [-1, 1] ...
                                                        );

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end