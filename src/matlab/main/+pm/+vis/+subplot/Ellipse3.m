classdef Ellipse < pm.matlab.Handle
    %
    %   This is the Ellipse class for generating
    %   instances of 3-dimensional Ellipse plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       gramian
    %
    %           The MATLAB (function handle returning a) scalar, matrix, or 3D
    %           rectangle of shape ``[ndim, ndim]`` or ``[ndim, ndim, nell]``,
    %           each subset ``gramian(1:ndim, 1:ndim, iell)`` of which represents
    %           the ``iell`` Gramian matrix of a ndim-dimensional ellipsoid to visualize
    %           where ``ndim`` refers to the number of dimensions of the hyper-ellipsoids
    %           represented by the input arguments and ``nell`` is the number of such ellipsoids.
    %
    %           If the specified ``gramian`` is of shape ``[ndim, ndim]``,
    %           it will be transformed to ``repmat(gramian, 1, nell)``
    %           to take the shape ``[ndim, ndim, nell]`` before usage.
    %
    %           If the specified ``gramian`` is a scalar,
    %           it will be transformed to ``repmat(gramian, ndim, nell)``
    %           to take the shape ``[ndim, ndim, nell]`` before usage.
    %
    %           \note
    %
    %               While it is possible to pass a data directly to this class,
    %               it is highly recommended to pass a function handle that returns
    %               such data when called, allowing the visualization data to be
    %               dynamically updated when needed.
    %
    %       center
    %
    %           The input MATLAB (function handle returning a)
    %           vector of doubles of size ``ndim`` or matrix of shape ``[ndim, nell]``
    %           containing the coordinates of the centers of the ellipsoids to display
    %           where ``ndim`` refers to the number of dimensions of the hyper-ellipsoids
    %           represented by the input arguments and ``nell`` is the number of such ellipsoids.
    %
    %           If the specified ``center`` is of shape ``[ndim, 1]``
    %           it will be transformed to ``repmat(center, 1, nell)``
    %           to take the shape ``[ndim, nell]`` before usage.
    %
    %           If the specified ``center`` is of shape ``[1, nell]``
    %           it will be transformed to ``repmat(center, ndim, 1)``
    %           to take the shape ``[ndim, nell]`` before usage.
    %
    %           \note
    %
    %               While it is possible to pass data directly to this class,
    %               it is highly recommended to pass a function handle that returns
    %               such data when called, allowing the visualization data to be
    %               dynamically updated when needed.
    %
    %       zval
    %
    %           The input MATLAB (function handle returning a)
    %           scalar, or vector of size ``npnt``, or matrix of shape ``[npnt, nell]``
    %           of doubles containing the z-axis coordinates of each of the ellipsoid projections
    %           along a given pair of dimensions specified via the ``dims`` attribute of parent object.
    %           where ``nell`` is the number of Gramian matrices in the input ``gramian`` and ``npnt``
    %           is the number of points used in visualizing each of the ellipsoid projections.
    %           (**optional**. The default is the ``ellindex`` attributes of the object.)
    %
    %           If the specified ``center`` is of shape ``[ndim, 1]``
    %           it will be transformed to ``repmat(center, 1, nell)``
    %           to take the shape ``[ndim, nell]`` before usage.
    %
    %           If the specified ``center`` is of shape ``[1, nell]``
    %           it will be transformed to ``repmat(center, ndim, 1)``
    %           to take the shape ``[ndim, nell]`` before usage.
    %
    %           If the specified ``zval`` is a scalar
    %           it will be transformed to ``zval * ones(zval, 1, nell)``
    %           to take the shape ``[npnt, nell]`` before usage.
    %
    %           \note
    %
    %               While it is possible to pass data directly to this class,
    %               it is highly recommended to pass a function handle that returns
    %               such data when called, allowing the visualization data to be
    %               dynamically updated when needed.
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``make()`` method.
    %
    %   Attributes
    %   ----------
    %
    %       See below and the documentation of the attributes
    %       of the parent class ``pm.vis.subplot.LineScatter3``.
    %
    %   Returns
    %   -------
    %
    %       An object of class ``pm.vis.subplot.Ellipse``.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.Ellipse(gramian, center);
    %       p = pm.vis.subplot.Ellipse(gramian, center, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties (Access = public)
        %
        %   npnt
        %
        %       The scalar MATLAB whole number representing the number of
        %       points with which the 2D ellipsoid boundaries are delineated.
        %       The default value is ``size(zval, 1)`` or otherwise, ``50``.
        %
        npnt = [];
        %
        %   dims
        %
        %       The vector of size ``2`` of MATLAB whole numbers representing the dimensions
        %       of the input ``gramian`` and ``center`` data to use for 2D ellipsoid visualization.
        %       The default is ``[1, 2]``.
        %
        dims = [];
        %
        %   ellindex
        %
        %       The vector of MATLAB whole numbers representing the ellindex of the 2D
        %       ellipsoids to display, represented by the input ``gramian`` and ``center``.
        %       The specified ellindex will serve as the visualization values on the z-axis.
        %       The default is ``ellindex = pm.array.logrange(start = 1, stop = nell, count = 100);``,
        %       where ``nell`` represents the number of ellipsoids to visualize.
        %
        ellindex = [];
        %
        %   names
        %
        %       The vector of MATLAB strings containing the names of dimensions
        %       of the space within which the Gramian matrices are defined.
        %       The default is ``"Dimension i"`` where ``i`` is replaced
        %       by the ID of the corresponding axis.
        %
        names = [];
        %
        %   gramian
        %
        %       A scalar object of class ``pm.data.DataRef``
        %       containing the user-specified Gramian data to visualize.
        %
        gramian = [];
        %
        %   center
        %
        %       A scalar object of class ``pm.data.DataRef``
        %       containing the user-specified center data to visualize.
        %
        center = [];
        %
        %   zval
        %
        %       A scalar object of class ``pm.data.DataRef``
        %       containing the user-specified z-axis data to visualize.
        %
        zval = [];
        %
        %       cval
        %
        %           A scalar object of class ``pm.data.DataRef`` whose value is a string,
        %           or scalar, or vector of size ``npnt``, or matrix of shape ``[npnt, nell]``
        %           of doubles, each subset ``(:, iell)`` containing the color of each of the ellipsoid projections
        %           along a given pair of dimensions specified via the ``dims`` attribute of parent object.
        %           where ``nell`` is the number of Gramian matrices in the input ``gramian`` and ``npnt``
        %           is the number of points used in visualizing each of the ellipsoid projections.
        %           (**optional**. The default is the ``ellindex`` attributes of the object.)
        %
        %           \note
        %
        %               While it is possible to pass data directly to this class,
        %               it is highly recommended to pass a function handle that returns
        %               such data when called, allowing the visualization data to be
        %               dynamically updated when needed.
        %
        cval = [];
    end

    properties (Access = public, Hidden)
        %
        %   ndim
        %
        %       A scalar MATLAB whole number representing the total number of the
        %       dimensions of ellipsoids identified based on the input that can be visualized.
        %       It is set to ``max(2, size(gramian, 2), size(center, 1))``.
        %
        ndim = [];
        %
        %   nell
        %
        %       A scalar MATLAB whole number representing the total number of
        %       ellipsoids identified based on the input that can be visualized.
        %       It is set to ``max(20, size(gramian, 3), size(center, 2), size(zval, 2))``.
        %
        nell = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Ellipse(gramian, center, zval, varargin)

            %%%% Define the missing optional values.

            if  nargin < 3
                zval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end

            %%%% Convert all data to dataframe reference temporarily.

            zval_ = pm.data.DataRef(zval, "VariableNames", "val");
            center_ = pm.data.DataRef(center, "VariableNames", "val");
            gramian_ = pm.data.DataRef(gramian, "VariableNames", "val");

            gramian_ = gramian_.copy();
            center_ = center_.copy();
            zval_ = zval_.copy();

            ndim = max(2, size(gramian_, 2), size(center_, 1));
            nell = max(20, size(gramian_, 3), size(center_, 2), size(zval_, 2));

            %%%% Call the parent constructor.

            varargin = {"ndim", ndim, "nell", nell, varargin{:}};
            self = self@pm.vis.subplot.LineScatter3([], varargin{:});

            %%%% Assert shape consistencies.

            asserted = true;
            asserted = asserted && ndim == size(gramian_, 2) || size(gramian_, 2) == 0;
            asserted = asserted && ndim == size(center_, 1) || size(center_, 1) == 0;
            asserted = asserted && ndim == size(zval_, 1) || size(zval_, 1) == 0;
            asserted = asserted && nell == size(gramian_, 3) || size(gramian_, 3) == 0;
            asserted = asserted && nell == size(center_, 2) || size(center_, 2) == 0;
            asserted = asserted && nell == size(zval_, 2) || size(zval_, 2) == 0;

            if ~asserted
                help("pm.vis.subplot.Ellipse")
                disp("size(gramian)")
                disp( size(gramian_) )
                disp("size(center)")
                disp( size(center_) )
                disp("size(zval)")
                disp( size(zval_) )
                disp("[ndim, nell]")
                disp( [ndim, nell] )
                error   ( newline ...
                        + "The shapes of the specified ``gramian``, ``center``, and ``zval`` are incompatible." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%% Redefine the empty input values.

            if  isempty(gramian_)
                gramian = eye(ndim, ndim, nell);
            end
            if  isempty(center_)
                center = zeros(ndim, nell);
            end
            if  isempty(zval_)
                zval = repmat(1 : nell, npnt, 1);
            end

            self.gramian = pm.data.DataRef(gramian);
            self.center = pm.data.DataRef(center);
            self.zval = pm.data.DataRef(zval);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
            %
            %   Generate a plot from the selected dimensions
            %   of the input ellipsoid data to the object constructor.
            %
            %   \note
            %
            %       This method is pure, except for the changes in ``fout`` component.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes, before calling the ``make()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.subplot.Subplot.make(varargin);
            %
            %   Example
            %   -------
            %
            %       dfref = rand(1000, 3);
            %       p = pm.vis.subplot.Subplot("scatter", dfref);
            %       p.make("colx", 1, "coly", 2, "colc", 3)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            premake@pm.vis.LineScatter3(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%% Set the data.

            zval = self.zval.copy();
            center = self.center.copy();
            gramian = self.gramian.copy();

           %%%% Check the consistency of ``self.dims``.

            if  length(self.dims) ~= 2
                help("pm.vis.subplot.Ellipse")
                disp("self.dims")
                disp( self.dims )
                disp("size(self.dims)")
                disp( size(self.dims) )
                error   ( newline ...
                        + "The size of the specified ``dims`` must be ``2``." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end

           %%%% Check the consistency of ``self.names``.

            if  length(self.names) ~= self.ndim
                help("pm.vis.subplot.Ellipse")
                disp("string(self.names)")
                disp( string(self.names) )
                disp("length(string(self.names))")
                disp( length(string(self.names)) )
                error   ( newline ...
                        + "The component ``names`` must be a vector of size ``ndim`` of strings." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%% Set up dimension df.

            if ~isempty(self.npnt)
                npnt = self.npnt;
            else
                if  1 < size(zval, 1)
                    npnt = size(zval, 1);
                else
                    npnt = 50;
                end
            end

            zvalScalar = numel(zval) == 1;
            bcrd = zeros(npnt, 3);
            for iell = 1 : size(self.ellindex)
                icol = (iell - 1) * 3 + 1;
                bcrd(:, icol : icol + 1) = pm.geom.ell2.getBorder(gramian(:, :, min(iell, size(gramian, 3))), center(:, min(iell, size(center, 2))), [], npnt);
                if  zvalScalar
                    bcrd(:, icol + 2) = zval;
                else
                    bcrd(:, icol + 2) = zval(:, min(iell, size(zval, 2)));
                end
            end

            self.df = pm.data.DataFrame()

            make@pm.vis.LineScatter3(self);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function premake(self, varargin)
            %
            %   Prepare the subplot for visualization.
            %
            %   \warning
            %
            %       This method has side-effects by manipulating
            %       the existing attributes of the parent object.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes, before calling the ``premake()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.subplot.Ellipse.premake(varargin);
            %
            %   Example
            %   -------
            %
            %       p = pm.vis.subplot.Ellipse();
            %       p.premake("npnt", 50, "dims", [1, 2])
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            premake@pm.vis.LineScatter3(self, varargin{:});
            if  isempty(self.dims)
                self.dims = [1, 2];
            end
            %if  isempty(self.npnt)
            %    if  1 < size(self.zval, 1)
            %        self.npnt = size(self.zval, 1);
            %    else
            %        self.npnt = [];
            %    end
            %end
            if ~isempty(self.names)
                names = string(self.names);
                names = names(:);
            else
                names = strings(self.ndim, 1);
                for idim = 1, self.ndim
                    self.names(idim) = "Dimension " + string(idim);
                end
            end
            if  isempty(self.ellindex)
                self.ellindex = pm.array.logrange(1, self.nell);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetint(self, varargin)
            self.dims = [];
            self.npnt = [];
            self.names = [];
            self.ellindex = [];
            resetint@pm.vis.LineScatter3(self, varargin{:});
            self.cmap.map = "winter";
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end