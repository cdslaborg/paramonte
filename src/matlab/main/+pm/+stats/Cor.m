%>  \brief
%>  This is the base class for generating objects with methods
%>  and storage components for computing and storing the
%>  correlation matrix of an input data.<br>
%>
%>  \details
%>  This is convenience class for easy computation
%>  of correlation and its storage all in one place.<br>
%>  The primary advantage of this class over the MATLAB
%>  intrinsic functions is in the ability of this class
%>  to compute the result for input dataframe table and
%>  return the results always in MATLAB ``table`` format.<br>
%>
%>  \note
%>  See the documentation of the class constructor below.<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:16 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Cor < pm.matlab.Handle

    properties(Access = public)
        %>
        %>  ``method``
        %>
        %>  The scalar MATLAB string containing the
        %>  method of computing the correlation matrix.<br>
        %>  It can be either:<br>
        %>  <ol>
        %>      <li>    ``"pearson"``   : for computing the Pearson correlation matrix of the input data.
        %>      <li>    ``"kendall"``   : for computing the kendall rank correlation matrix of the input data.
        %>      <li>    ``"spearman"``  : for computing the Spearman rank correlation matrix of the input data.
        %>  </ol>
        %>
        method = "pearson";
        %>
        %>  ``val``
        %>
        %>  The MATLAB table of rank ``2`` serving as a
        %>  convenient storage component for the correlation matrix.<br>
        %>  This component is automatically populated at the time of
        %>  constructing an object of class [pm.stats.Cor](@ref Cor).<br>
        %>  It must be populated manually at all other times.<br>
        %>
        val = [];
    end

    methods(Access=public)

        %>  \brief
        %>  Return an object of class ``Cor``.<br>
        %>
        %>  \details
        %>  This is the constructor of the ``Cor`` class.<br>
        %>
        %>  \param[in]  df      :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose correlation matrix must be computed.<br>
        %>                          (**optional**. If missing, the correlation matrix will not be computed.)
        %>
        %>  \param[in]  method  :   The input scalar MATLAB string that can be either:<br>
        %>                          <ol>
        %>                              <li>    ``"pearson"``   : for computing the Pearson correlation matrix of the input data.<br>
        %>                              <li>    ``"kendall"``   : for computing the kendall rank correlation matrix of the input data.<br>
        %>                              <li>    ``"spearman"``  : for computing the Spearman rank correlation matrix of the input data.<br>
        %>                          </ol>
        %>                          (**optional**, default = ``"pearson"``)
        %>
        %>  \return
        %>  ``self``            :   The output object of class [pm.stats.Cor](@ref Cor).
        %>
        %>  \interface{Cor}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor()
        %>      mat = pm.stats.Cor(df)
        %>      mat = pm.stats.Cor(df, method)
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
        function self = Cor(df, method)
            if nargin == 2
                self.method = method;
            end
            if  0 < nargin
                self.val = self.get(df, self.method);
            end
        end

        %>  \brief
        %>  Return the correlation matrix of the input data.<br>
        %>
        %>  \details
        %>  This is a dynamic method of the ``Cor`` class.
        %>  This method automatically stores any input information
        %>  in the corresponding components of the parent object.<br>
        %>  However, any components of the parent object
        %>  corresponding to the output of this method
        %>  must be set explicitly manually.<br>
        %>
        %>  \param[in]  df      :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose correlation matrix must be computed.<br>
        %>
        %>  \param[in]  method  :   The input scalar MATLAB string that can be either:<br>
        %>                          <ol>
        %>                              <li>    ``"pearson"``   : for computing the Pearson correlation matrix of the input data.<br>
        %>                              <li>    ``"kendall"``   : for computing the kendall rank correlation matrix of the input data.<br>
        %>                              <li>    ``"spearman"``  : for computing the Spearman rank correlation matrix of the input data.<br>
        %>                          </ol>
        %>                          (**optional**, default = [pm.stats.Cor.method](@ref Cor::method))
        %>
        %>  \return
        %>  ``val``             :   The output MATLAB ``table`` containing the correlation matrix.<br>
        %>
        %>  \interface{get}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor()
        %>      mat.val = mat.get(df)
        %>      mat.val = mat.get(df, method)
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
        function val = get(self, df, method)
            if  nargin < 2
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
            catch me
                val = NaN(size(data, 2), size(data, 2));
                warning ( newline ...
                        + string(me.identifier) + " : " + string(me.message) + newline ...
                        + "skipping the correlation matrix computation..." + newline ...
                        );
                return;
            end
            if isa(df, "table")
                val.Properties.VariableNames = df.Properties.VariableNames;
            end
            val.Properties.RowNames = val.Properties.VariableNames;
        end

    end

end