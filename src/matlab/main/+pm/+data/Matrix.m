classdef Matrix < pm.matlab.Handle
    %
    %   This is the base class for generating objects with methods
    %   and storage components for storing and visualizing square
    %   positive-definite matrices
    %
    %   \note
    %
    %       This is convenience class for easy storage and visualization
    %       of square positive-definite matrices all in one place.
    %       The primary advantage of this class over the MATLAB
    %       intrinsic functions is in the ability of this class
    %       to visualize the result for input dataframe table and
    %       retain the results always in MATLAB ``table`` format.
    %
    %   Parameters
    %   ----------
    %
    %       See the documentation of the class constructor below.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the class constructor below.
    %
    %   Returns
    %   -------
    %
    %       See the documentation of the class constructor below.
    %
    %   Interface
    %   ---------
    %
    %       matrix = pm.data.Matrix()
    %       matrix = pm.data.Matrix(df)
    %       matrix = pm.data.Matrix(df, method)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties(Access = public)
        %
        %   val
        %
        %       The MATLAB table of rank ``2`` serving as a
        %       convenient storage component for the square positive-definite matrix.
        %       This component is automatically populated at the time of
        %       constructing an object of class ``pm.data.Matrix``.
        %       It can be populated manually at all other times.
        %
        val = [];
    end

    methods(Access = public)

        function self = matrix(val)
            %
            %   Return an object of class ``matrix``.
            %
            %   This is the constructor of the ``Matrix`` class.
            %
            %   Parameters
            %   ----------
            %
            %       val
            %
            %           The input MATLAB matrix or table of rank ``2``
            %           containing the square positive matrix with ``ndim`` rows and ``ndim columns.
            %           (**optional**. If missing, the ``val`` component of the output object will be empty.)
            %
            %   Returns
            %   -------
            %
            %       self
            %
            %           The output object of class ``pm.data.Matrix``.
            %
            %   Interface
            %   ---------
            %
            %       matrix = pm.data.Matrix()
            %       matrix = pm.data.Matrix(df)
            %       matrix = pm.data.Matrix(df, method)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if  0 < nargin
                if size(val, 1) ~= size(val, 2)
                    error   ( newline ...
                            + "The input ``val`` must be a square positive-definite matrix or table." + newline ...
                            + newline ...
                            );
                end
                self.val = val;
            end
        end

        function val = get(self, df, method)
            %
            %   Return the covariance matrix of the input data.
            %
            %   This is a dynamic method of the ``Matrix`` class.
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
            %           observations whose covariance matrix must be computed.
            %
            %       method
            %
            %           The input scalar MATLAB string that can be either:
            %
            %               "pearson"   : for computing the Pearson covariance matrix of the input data.
            %               "spearman"  : for computing the Spearman rank covariance matrix of the input data.
            %
            %           (**optional**, default = ``"pearson"``)
            %
            %   Returns
            %   -------
            %
            %       val
            %
            %           The output MATLAB ``table`` containing the covariance matrix.
            %
            %   Interface
            %   ---------
            %
            %       matrix = pm.data.Matrix()
            %       matrix.val = matrix.get(df)
            %       matrix.val = matrix.get(df, method)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if nargin < 2
                help("pm.data.Matrix");
                error   ( newline ...
                        + "The input ``df`` argument is required for computing the covariance matrix." + newline ...
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
            if strcmpi(self.method, "spearman")
                try
                    data = tiedrank(data);
                catch ME
                    val = NaN(size(data, 2), size(data, 2));
                    warning ( newline ...
                            + string(ME.message) + newline ...
                            + newline ...
                            + "skipping the covariance matrix computation..." + newline ...
                            );
                    return;
                end
            end
            try
                val = array2table(cov(data));
            catch ME
                val = NaN(size(data, 2), size(data, 2));
                warning ( newline ...
                        + string(ME.message) + newline ...
                        + newline ...
                        + "skipping the covariance matrix computation..." + newline ...
                        );
                return;
            end
            if  isa(df, "table")
                val.Properties.VariableNames = df.Properties.VariableNames;
            end
        end

    end

end