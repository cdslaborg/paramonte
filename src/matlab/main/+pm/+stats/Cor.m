%>  \brief
%>  This is the base class for generating objects with methods
%>  and storage components for computing and storing the
%>  correlation matrix of an input data.
%>
%>  \note
%>  This is convenience class for easy computation
%>  of correlation and its storage all in one place.<br>
%>  The primary advantage of this class over the MATLAB
%>  intrinsic functions is in the ability of this class
%>  to compute the result for input dataframe table and
%>  return the results always in MATLAB ``table`` format.<br>
%>  See the documentation of the class constructor below.
%>
%>  \note
%>  See the documentation of the class constructor below.
%>
%>  \return
%>  See the documentation of the class constructor below.
%>
%>  \interface{Cor}
%>  \code{.m}
%>
%>      mat = pm.stats.Cor()
%>      mat = pm.stats.Cor(df)
%>      mat = pm.stats.Cor(df, method)
%>
%>  \endcode
%>  \final{Cor}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:16 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Cor < pm.matlab.Handle

    properties(Access = public)
        %>  
        %>  \param  method  :   The scalar MATLAB string containing the
        %>                      method of computing the correlation matrix.
        %>                      It can be either:<br>
        %>                      "pearson"   : for computing the Pearson correlation matrix of the input data.
        %>                      "kendall"   : for computing the kendall rank correlation matrix of the input data.
        %>                      "spearman"  : for computing the Spearman rank correlation matrix of the input data.
        %>  
        %>  
        method = "pearson";
        %>
        %>  \param  val     :   The MATLAB table of rank ``2`` serving as a
        %>                      convenient storage component for the correlation matrix.<br>
        %>                      This component is automatically populated at the time of
        %>                      constructing an object of class ``pm.stats.Cor``.<br>
        %>                      It must be populated manually at all other times.
        %>
        val = [];
    end

    methods(Access=public)

        %>  \brief
        %>  Return an object of class ``Cor``.<br>
        %>  This is the constructor of the ``Cor`` class.
        %>
        %>  \param[in]  df      :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose correlation matrix must be computed.<br>
        %>                          (**optional**. If missing, the correlation matrix will not be computed.)
        %>  
        %>  \param[in]  method  :   The input scalar MATLAB string that can be either:<br>
        %>                          "pearson"   : for computing the Pearson correlation matrix of the input data.<br>
        %>                          "kendall"   : for computing the kendall rank correlation matrix of the input data.<br>
        %>                          "spearman"  : for computing the Spearman rank correlation matrix of the input data.<br>
        %>                          (**optional**, default = ``"pearson"``)
        %>
        %>  \return
        %>  `self`              :   The output object of class ``pm.stats.Cor``.
        %>
        %>  \interface{Cor}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor()
        %>      mat = pm.stats.Cor(df)
        %>      mat = pm.stats.Cor(df, method)
        %>
        %>  \endcode
        %>  \final{Cor}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:22 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Cor(df, method)
            if nargin == 2
                self.method = method;
            end
            if 0 < nargin
                self.val = self.get(df, self.method);
            end
        end

        %>  \brief
        %>  Return the correlation matrix of the input data.
        %>
        %>  \details
        %>  This is a dynamic method of the ``Cor`` class.
        %>  This method automatically stores any input information
        %>  in the corresponding components of the parent object.
        %>  However, any components of the parent object
        %>  corresponding to the output of this method
        %>  must be set explicitly manually.
        %>
        %>  \param[in]  df      :   The input MATLAB matrix or table of rank ``2``
        %>                          containing the data as ``ncol`` columns of ``nrow``
        %>                          observations whose correlation matrix must be computed.
        %>  
        %>  \param[in]  method  :   The input scalar MATLAB string that can be either:<br>
        %>                          "pearson"   : for computing the Pearson correlation matrix of the input data.
        %>                          "kendall"   : for computing the kendall rank correlation matrix of the input data.
        %>                          "spearman"  : for computing the Spearman rank correlation matrix of the input data.<br>
        %>                          (**optional**, default = ``pm.stats.Cor.method``)
        %>
        %>  \return
        %>  `val`               :   The output MATLAB ``table`` containing the correlation matrix.
        %>
        %>  \interface{get}
        %>  \code{.m}
        %>
        %>      mat = pm.stats.Cor()
        %>      mat.val = mat.get(df)
        %>      mat.val = mat.get(df, method)
        %>
        %>  \endcode
        %>  \final{get}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 4:24 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
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
        end

    end

end