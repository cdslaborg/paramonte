%>  \brief
%>  This is the base class for generating objects with methods and
%>  storage components for computing and storing and visualizing
%>  the covariance matrix of an input data.<br>
%>
%>  \details
%>  This is merely a convenience class for easy computation of covariance and its storage all in one place.<br>
%>  The primary advantage of this class over the MATLAB intrinsic functions is in the ability of this class to
%>  compute the result for input dataframe table and return the results always in MATLAB ``table`` format.<br>
%>
%>  \note
%>  See the documentation of the class constructor [pm.stats.Cov](@ref Cov::Cov) below.<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:25 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
classdef Cov < pm.matlab.Handle

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
        %>  ``method``
        %>
        %>  The scalar MATLAB string containing the
        %>  method of computing the covariance matrix.<br>
        %>  It can be either:<br>
        %>  <ol>
        %>      <li>    ``"pearson"``   : for computing the Pearson covariance matrix of the input data.<br>
        %>      <li>    ``"spearman"``  : for computing the Spearman rank covariance matrix of the input data.<br>
        %>  </ol>
        %>
        method = "pearson";
        %>
        %>  ``val``
        %>
        %>  The MATLAB ``table`` of rank ``2`` serving as a
        %>  convenient storage component for the covariance matrix.<br>
        %>  This component is automatically populated at the time of
        %>  constructing an object of class [pm.stats.Cov](@ref Cov).<br>
        %>  It must be populated manually at all other times.<br>
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Return an object of class [pm.stats.Cov](@ref Cov).<br>
        %>
        %>  \details
        %>  This is the constructor of the [pm.stats.Cov](@ref Cov) class.<br>
        %>
        %>  \param[in]  dfref   :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose covariance matrix must be computed.<br>
        %>                          Ideally, the user would want to pass a reference
        %>                          to a dataframe (e.g., as a function handle ``@()df``)
        %>                          so that the data remains dynamically up-to-date.<br>
        %>                          (**optional**. If missing or empty, the covariance matrix will not be computed.)
        %>  \param[in]  method  :   The input scalar MATLAB string that can be either:<br>
        %>                          <ol>
        %>                              <li>    ``"pearson"``   : for computing the Pearson covariance matrix of the input data.<br>
        %>                              <li>    ``"spearman"``  : for computing the Spearman rank covariance matrix of the input data.<br>
        %>                          </ol>
        %>                          (**optional**, default = ``"pearson"``)
        %>
        %>  \return
        %>  ``self``            :   The output object of class [pm.stats.Cov](@ref Cov).<br>
        %>
        %>  \interface{Cov}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cov([])
        %>      mat = pm.stats.Cov([], [])
        %>      mat = pm.stats.Cov([], method)
        %>
        %>      mat = pm.stats.Cov(dfref)
        %>      mat = pm.stats.Cov(dfref, [])
        %>      mat = pm.stats.Cov(dfref, method)
        %>
        %>  \endcode
        %>
        %>  \example{Cov}
        %>  \include{lineno} example/stats/Cov/main.m
        %>  \output{Cov}
        %>  \include{lineno} example/stats/Cov/main.out.m
        %>  \vis{Cov}
        %>  \image html example/stats/Cov/Cov.unifrnd.png width=700
        %>
        %>  \final{Cov}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:29 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Cov(dfref, method)

            if  1 < nargin
                self.method = method;
            end

            if  0 < nargin
                self.val = self.get(dfref, self.method);
            end

            self.setvis();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Return the covariance matrix of the input data.<br>
        %>
        %>  \details
        %>  This is a dynamic method of the [pm.stats.Cov](@ref Cov) class.<br>
        %>  This method automatically stores any input information in
        %>  the corresponding components of the parent object.<br>
        %>  However, any components of the parent object
        %>  corresponding to the output of this method
        %>  must be set explicitly manually.<br>
        %>
        %>  \param[inout]   self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]      dfref   :   The input (reference of function handle returning a) MATLAB matrix or table of rank ``2``
        %>                              containing the ``ncol`` columns of ``nrow`` data whose covariance matrix must be computed.<br>
        %>                              Ideally, the user would want to pass a reference to a dataframe (e.g., as a function handle ``@()df``)
        %>                              so that the data remains dynamically up-to-date.<br>
        %>                              (**optional**. If missing, the contents of the corresponding internal component of the parent object will be used.)
        %>  \param[in]      method  :   The input scalar MATLAB string that can be either:<br>
        %>                              <ol>
        %>                                  <li>    ``"pearson"``   : for computing the Pearson covariance matrix of the input data.<br>
        %>                                  <li>    ``"spearman"``  : for computing the Spearman rank covariance matrix of the input data.<br>
        %>                              </ol>
        %>                              (**optional**, default = ``"pearson"``)
        %>
        %>  \return
        %>  ``val``                 :   The output MATLAB ``table`` containing the covariance matrix.<br>
        %>
        %>  \interface{get}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cov(dfref, method)
        %>
        %>      mat.val = mat.get()
        %>      mat.val = mat.get([])
        %>      mat.val = mat.get([], [])
        %>      mat.val = mat.get([], method)
        %>
        %>      mat.val = mat.get(dfref)
        %>      mat.val = mat.get(dfref, [])
        %>      mat.val = mat.get(dfref, method)
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the documentation of the class constructor
        %>  [pm.stats.Cor](@ref Cor::Cor) for example usage.<br>
        %>
        %>  \final{get}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:31 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function val = get(self, dfref, method)

            if  nargin < 3
                method = "pearson";
            end
            self.method = method;

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
                help("pm.stats.Cov");
                error   ( newline ...
                        + "A non-mepty ``dfref`` attribute or input argument is required for computing the covariance matrix." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Convert data to rank if necessary.
            %%%%

            if  strcmpi(self.method, "spearman")
                try
                    data = tiedrank(data);
                catch me
                    val = NaN(size(data, 2), size(data, 2));
                    warning ( newline ...
                            + string(me.identifier) + " : " + string(me.message) + newline ...
                            + "skipping the covariance matrix computation..." + newline ...
                            + newline ...
                            );
                    return;
                end
            elseif ~strcmpi(self.method, "pearson")
                help("pm.stats.Cov");
                disp("self.method");
                disp( self.method );
                error   ( newline ...
                        + "Invalid value specified for the covariance computation ``method``." + newline ...
                        + newline ...
                        );
            end

            try
                val = array2table(cov(data));
            catch me
                val = NaN(size(data, 2), size(data, 2));
                warning ( newline ...
                        + string(me.identifier) + " : " + string(me.message) + newline ...
                        + "skipping the covariance matrix computation..." + newline ...
                        + newline ...
                        );
                return;
            end

            %%%%
            %%%% Set the matrix column and row names.
            %%%%

            if  isa(dfcopy, "table")
                val.Properties.VariableNames = dfcopy.Properties.VariableNames;
            end
            val.Properties.RowNames = val.Properties.VariableNames;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Set up the visualization tools of the covariance matrix.<br>
        %>
        %>  \details
        %>  This is a dynamic ``Hidden`` method of the [pm.stats.Cov](@ref Cov) class.<br>
        %>  This method is inaccessible to the end users of the ParaMonte MATLAB library.<br>
        %>
        %>  \param[inout]   self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]      val     :   The input (reference of function handle returning a) MATLAB matrix or table of rank ``2``
        %>                              containing the computed covariance matrix to be visualized.<br>
        %>                              Ideally, the user would want to pass a reference to a dataframe (e.g., as a function handle ``@()df``)
        %>                              so that the data remains dynamically up-to-date.<br>
        %>                              (**optional**. If missing, the contents of the corresponding ``var`` attribute of the parent object will be used.)
        %>
        %>  \interface{setvis}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cov(dfref, method)
        %>      mat.setvis(); % This method is automatically called within the object constructor.
        %>
        %>  \endcode
        %>
        %>  \final{setvis}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function setvis(self, val)

            if  1 < nargin
                self.val = val;
            end

            self.vis = struct();
            self.vis.heatmap = pm.vis.PlotHeatmap(@()self.val);
            self.vis.heatmap.subplot.title.titletext = "Covariance Matrix";
            self.vis.heatmap.subplot.title.enabled = true;
            self.vis.heatmap.subplot.precision = 2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end