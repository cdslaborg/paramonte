%>  \brief
%>  This is the base class for generating objects with methods
%>  and storage components for computing, storing, and visualizing
%>  the correlation matrix of an input data.<br>
%>
%>  \details
%>  This is convenience class for easy correlation computation
%>  and its storage and visualization all in one place.<br>
%>  The primary advantage of this class over the MATLAB
%>  intrinsic functions is in the ability of this class
%>  to compute the result for input dataframe table and
%>  return the results always in MATLAB ``table`` format.<br>
%>
%>  \note
%>  See the documentation of the class constructor below.<br>
%>
%>  \note
%>  See also the documentation of the superclass [pm.stats.Cov](@ref Cov).<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:16 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Cor < pm.stats.Cov

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Return an object of class [pm.stats.Cor](@ref Cor).<br>
        %>
        %>  \details
        %>  This is the constructor of the [pm.stats.Cor](@ref Cor) class.<br>
        %>
        %>  \param[in]  dfref   :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose correlation matrix must be computed.<br>
        %>                          Ideally, the user would want to pass a reference
        %>                          to a dataframe (e.g., as a function handle ``@()df``)
        %>                          so that the data remains dynamically up-to-date.<br>
        %>                          (**optional**. If missing or empty, the correlation matrix will not be computed.)
        %>  \param[in]  method  :   The input scalar MATLAB string that can be either:<br>
        %>                          <ol>
        %>                              <li>    ``"pearson"``   : for computing the Pearson correlation matrix of the input data.<br>
        %>                              <li>    ``"kendall"``   : for computing the Kendall correlation matrix of the input data.<br>
        %>                              <li>    ``"spearman"``  : for computing the Spearman rank correlation matrix of the input data.<br>
        %>                          </ol>
        %>                          (**optional**, default = ``"pearson"``)
        %>
        %>  \return
        %>  ``self``            :   The output object of class [pm.stats.Cor](@ref Cor).<br>
        %>
        %>  \interface{Cor}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor([])
        %>      mat = pm.stats.Cor([], [])
        %>      mat = pm.stats.Cor([], method)
        %>
        %>      mat = pm.stats.Cor(dfref)
        %>      mat = pm.stats.Cor(dfref, [])
        %>      mat = pm.stats.Cor(dfref, method)
        %>
        %>  \endcode
        %>
        %>  \example{Cor}
        %>  \include{lineno} example/stats/Cor/main.m
        %>  \output{Cor}
        %>  \include{lineno} example/stats/Cor/main.out.m
        %>  \vis{Cor}
        %>  \image html example/stats/Cor/Cor.unifrnd.png width=700
        %>
        %>  \final{Cor}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:22 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Cor(dfref, method)

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
        %>  Return the correlation matrix of the input data.<br>
        %>
        %>  \details
        %>  This is a dynamic method of the ``Cor`` class.<br>
        %>  This method automatically stores any input information
        %>  in the corresponding components of the parent object.<br>
        %>  However, any components of the parent object
        %>  corresponding to the output of this method
        %>  must be set explicitly manually.<br>
        %>
        %>  \param[inout]   self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]      dfref   :   The input (reference of function handle returning a) MATLAB matrix or table of rank ``2``
        %>                              containing the ``ncol`` columns of ``nrow`` data whose correlation matrix must be computed.<br>
        %>                              Ideally, the user would want to pass a reference to a dataframe (e.g., as a function handle ``@()df``)
        %>                              so that the data remains dynamically up-to-date.<br>
        %>                              (**optional**. If missing, the correlation matrix will not be computed.)
        %>  \param[in]      method  :   The input scalar MATLAB string that can be either:<br>
        %>                              <ol>
        %>                                  <li>    ``"pearson"``   : for computing the Pearson correlation matrix of the input data.<br>
        %>                                  <li>    ``"kendall"``   : for computing the kendall rank correlation matrix of the input data.<br>
        %>                                  <li>    ``"spearman"``  : for computing the Spearman rank correlation matrix of the input data.<br>
        %>                              </ol>
        %>                              (**optional**, default = ``"pearson"``)
        %>
        %>  \return
        %>  ``val``             :   The output MATLAB ``table`` containing the correlation matrix.<br>
        %>
        %>  \interface{get}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor(dfref, method)
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
        %>  \JoshuaOsborne, May 21 2024, 4:24 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function val = get(self, dfref, method)

            if  nargin < 3
                method = "pearson";
            end

            if  nargin < 2
                dfref = [];
            end

            self.method = method;
            if ~isempty(dfref)
                self.dfref = pm.container.DataFrame(dfref);
            end

            dfcopy = self.dfref.copy();

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
            %%%% Compute the correlation matrix.
            %%%%

            try
                val = array2table(corr(data, "type", self.method));
            catch me
                val = NaN(size(data, 2), size(data, 2));
                warning ( newline ...
                        + string(me.identifier) + " : " + string(me.message) + newline ...
                        + "skipping the correlation matrix computation..." + newline ...
                        );
                return;
            end

            %%%%
            %%%% Set the matrix column and row names.
            %%%%

            if isa(dfcopy, "table")
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
        %>  Set up the visualization tools of the correlation matrix.<br>
        %>
        %>  \details
        %>  This is a dynamic ``Hidden`` method of the [pm.stats.Cor](@ref Cor) class.<br>
        %>  This method is inaccessible to the end users of the ParaMonte MATLAB library.<br>
        %>
        %>  \param[inout]   self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]      val     :   The input (reference of function handle returning a) MATLAB matrix or table of rank ``2``
        %>                              containing the computed correlation matrix to be visualized.<br>
        %>                              Ideally, the user would want to pass a reference to a dataframe (e.g., as a function handle ``@()df``)
        %>                              so that the data remains dynamically up-to-date.<br>
        %>                              (**optional**. If missing, the contents of the corresponding ``var`` attribute of the parent object will be used.)
        %>
        %>  \interface{setvis}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor(dfref, method)
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

            setvis@pm.stats.Cov(self);
            self.vis.heatmap.subplot.title.titletext = "Correlation Matrix";

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end