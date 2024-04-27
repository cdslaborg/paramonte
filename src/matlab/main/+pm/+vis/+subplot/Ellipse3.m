classdef Ellipse3 < pm.vis.subplot.LineScatter3
    %
    %   This is the Ellipse3 class for generating
    %   instances of 3-dimensional Ellipse3 plots
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
    %       An object of class ``pm.vis.subplot.Ellipse3``.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.Ellipse3();
    %       p = pm.vis.subplot.Ellipse3(gramian);
    %       p = pm.vis.subplot.Ellipse3(gramian, center);
    %       p = pm.vis.subplot.Ellipse3(gramian, center, val);
    %       p = pm.vis.subplot.Ellipse3(gramian, center, val, varargin);
    %
    %   Example
    %   -------
    %
    %       p = pm.vis.subplot.Ellipse3();
    %       p.make("npnt", 200);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties(Access = public)
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
        %   npnt
        %
        %       The scalar MATLAB whole number representing the number of
        %       points with which the 2D ellipsoid boundaries are delineated.
        %       The default value is ``size(zval, 1)`` or otherwise, ``100``.
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
        %       The default is ``ellindex = pm.array.logrange(start = 1, stop = nell, count = 150);``,
        %       where ``nell`` represents the number of ellipsoids to visualize.
        %
        ellindex = [];
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
        %   cval
        %
        %       A scalar object of class ``pm.data.DataRef``
        %       containing the user-specified color data in colormap.
        %       If empty or unspecified, it will be set to ``zval`` component.
        %
        cval = [];
    end

    properties(Access = public, Hidden)
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
        %
        %   workspace
        %
        %       A scalar MATLAB ``struct()`` containing the current visualization information.
        %       This information is updated with every call to the ``make()`` method.
        %       The contents of this component are temporary and must remain read-only.
        %       This component is solely provided for better insight into the internal
        %       workings of the ``make()`` method of the parent object and
        %       the resulting visualization.
        %
        workspace = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        function self = Ellipse3(gramian, center, zval, cval, varargin)

            %%%% Define the missing optional values as empty with the right rank.

            if  nargin < 4
                cval = zeros(0, 0);
            end
            if  nargin < 3
                zval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end

            varargin =  { "gramian", pm.data.DataRef(gramian) ...
                        , "center", pm.data.DataRef(center) ...
                        , "zval", pm.data.DataRef(zval) ...
                        , "cval", pm.data.DataRef(cval) ...
                        , varargin{:} ...
                        };
            self = self@pm.vis.subplot.LineScatter3([], varargin{:})

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        function reset(self, varargin)
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot and
            %   you want to clean up and go back to the default settings.
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
            %       pm.vis.subplot.Ellipse.reset() # reset the plot to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self.dims = [];
            self.npnt = [];
            self.names = [];
            self.ellindex = [];
            reset@pm.vis.subplot.LineScatter3(self, varargin{:});
            self.surface.lineWidth = [];
            self.scatter3.enabled = [];
            self.surface.enabled = [];
            self.plot3.enabled = [];
            self.colormap.map = []; %"winter"
            self.colormap.enabled = [];
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)
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
            %       pm.vis.subplot.Ellipse3.premake(varargin);
            %
            %   Example
            %   -------
            %
            %       p = pm.vis.subplot.Ellipse3();
            %       p.premake("npnt", 100, "dims", [1, 2])
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %

            %%%%
            %%%% Set the visualization types.
            %%%%

            if  isempty(self.colormap.map)
                self.colormap.map = "winter";
            end
            if  isempty(self.colormap.enabled)
                self.colormap.enabled = true;
            end

            if  isempty(self.plot3.enabled)
                self.plot3.enabled = false;
            end

            if  isempty(self.surface.enabled)
                self.surface.enabled = true;
            end

            if isempty(self.scatter3.enabled)
                self.scatter3.enabled = ~(self.plot3.enabled || self.surface.enabled);
            end

            %%%%
            %%%% Set the visualization specs.
            %%%%

            if  isempty(self.surface.lineWidth)
                self.surface.lineWidth = 1.5;
            end

            %%%%
            %%%% Set the visualization specs of the parent (and set the data components).
            %%%%

            premake@pm.vis.subplot.LineScatter3(self, varargin{:});

            %%%%
            %%%% Get the shape of input data.
            %%%%

            shape = struct();
            shape.gramian = size(self.gramian.copy(), 1, 2, 3);
            shape.center = size(self.center.copy(), 1, 2);
            shape.zval = size(self.zval.copy(), 1, 2);

            %%%%
            %%%% Infer ``ndim``.
            %%%%

            self.ndim = max(shape.gramian(1), shape.center(1));
            if  self.ndim == 0
                self.ndim = 2;
            end

            %%%%
            %%%% Infer ``nell``.
            %%%%

            self.nell = max([shape.gramian(3), shape.center(2), shape.zval(2)]);
            if  self.nell == 0
                self.nell = 20;
            end

            %%%%
            %%%% Infer ``npnt``.
            %%%%

            if  isempty(self.npnt)
                self.npnt = shape.zval(1);
                if  self.npnt == 0
                    self.npnt = 100;
                end
            end

            %%%%
            %%%% Assert data shape consistencies.
            %%%%

            asserted = true;
            asserted = asserted && self.ndim == shape.center(1) || shape.center(1) < 2;
            asserted = asserted && self.ndim == shape.gramian(1) || shape.gramian(1) < 2;
            asserted = asserted && self.ndim == shape.gramian(2) || shape.gramian(2) < 2;
            asserted = asserted && self.nell == shape.gramian(3) || shape.gramian(3) < 2;
            asserted = asserted && self.nell == shape.center(2) || shape.center(2) < 2;
            asserted = asserted && self.nell == shape.zval(2) || shape.zval(2) < 2;

            if ~asserted
                help("pm.vis.subplot.Ellipse3")
                disp("size(gramian)")
                disp( shape.gramian )
                disp("size(center)")
                disp( shape.center )
                disp("size(zval)")
                disp( shape.zval )
                disp("[ndim, nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``gramian``, ``center``, and ``zval`` are incompatible." + newline ...
                        + "For more information, see the documentation of the input arguments to the object constructor displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Set the default dims or verify the custom values to visualize.
            %%%%

            if  isempty(self.dims)
                self.dims = [1, 2];
            elseif  length(self.dims) ~= 2
                help("pm.vis.subplot.Ellipse3")
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

           %%%%
           %%%% Check the consistency of shape of ``self.names``.
           %%%%

            if ~isempty(self.names) && pm.array.len(self.names) ~= self.ndim
                help("pm.vis.subplot.Ellipse3")
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

            %%%%
            %%%% Set the data.
            %%%%

            self.workspace = struct();
            self.workspace.zval = self.zval.copy();
            self.workspace.center = self.center.copy();
            self.workspace.gramian = self.gramian.copy();

            %%%%
            %%%% Setup the gramian values.
            %%%%

            if  numel(self.workspace.gramian) == 0
                % The gramians of all ellipsoids are equal and identity matrix.
                self.workspace.gramian = repmat(eye(self.ndim, self.ndim), 1, 1, self.nell);
            elseif  numel(self.workspace.gramian) == 1
                % The gramians of all ellipsoids are equal and a multiple of the identity matrix.
                self.workspace.gramian = repmat(self.workspace.gramian * eye(self.ndim, self.ndim), 1, 1, self.nell);
            elseif  size(self.workspace.gramian, 1) == self.ndim && size(self.workspace.gramian, 2) == self.ndim && size(self.workspace.gramian, 3) == 1
                % All ellipsoids have the same user-prescribed self.workspace.gramian matrix.
                self.workspace.gramian = repmat(self.workspace.gramian(:, :, 1), 1, 1, self.nell);
            elseif  size(self.workspace.gramian, 1) == 1 && size(self.workspace.gramian, 2) == 1 && size(self.workspace.gramian, 3) == self.nell
                % Each ellipse has its unique diagonal self.workspace.gramian matrix whose diagonal values are all equal and prescribed by the user.
                temp = zeros(self.ndim, self.ndim, self.nell);
                for iell = 1 : self.nell
                    temp(:, :, iell) = self.workspace.gramian(1, 1, iell) * eye(self.ndim, self.ndim);
                end
                self.workspace.gramian = temp;
            elseif  size(self.workspace.gramian, 1) ~= self.ndim && size(self.workspace.gramian, 2) ~= self.ndim && size(self.workspace.gramian, 3) ~= self.nell
                help("pm.vis.subplot.Ellipse3")
                disp("size(gramian)")
                disp( size(self.workspace.gramian) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``gramian`` is incompatible with the" + newline ...
                        + "inferred number of ellipsoid dimensions ``ndim`` and the number of ellipsoids ``nell``." + newline ...
                        + "For more information, see the documentation of the ``gramian`` component displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Setup the center values.
            %%%%

            if  numel(self.workspace.center) == 0
                % The centers of all ellipsoids have the same coordinates and all coordinate axes have zero values.
                self.workspace.center = zeros(self.ndim, self.nell);
            elseif  numel(self.workspace.center) == 1
                % The centers of all ellipsoids have the same coordinates and all coordinate axes have equal non-zero values.
                self.workspace.center = self.workspace.center * ones(self.ndim, self.nell);
            elseif  size(self.workspace.center, 1) == self.ndim && size(self.workspace.center, 2) == 1
                % The centers of all ellipsoids have the same coordinates.
                self.workspace.center = repmat(self.workspace.center, 1, self.nell);
            elseif  size(self.workspace.center, 1) == 1 && size(self.workspace.center, 2) == self.nell
                % All coordinates of each ellipse self.workspace.center are the same.
                self.workspace.center = repmat(self.workspace.center, self.ndim, 1);
            elseif  size(self.workspace.center, 1) ~= self.ndim || size(self.workspace.center, 2) ~= self.nell
                help("pm.vis.subplot.Ellipse3")
                disp("size(self.workspace.center)")
                disp( size(self.workspace.center) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``center`` is incompatible with the" + newline ...
                        + "inferred number of ellipsoid dimensions ``ndim`` and the number of ellipsoids ``nell``." + newline ...
                        + "For more information, see the documentation of the ``center`` component displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Setup the zval values.
            %%%%

            if  numel(self.workspace.zval) == 0
                % 2D plot in 3D visualization at zero z value.
                self.workspace.zval = repmat(1 : self.nell, self.npnt, 1);
                if  isempty(self.axes.zscale)
                    self.axes.zscale = "log";
                end
            elseif  numel(self.workspace.zval) == 1
                % 2D plot in 3D visualization at a non-zero z value.
                self.workspace.zval = self.workspace.zval * ones(self.npnt, self.nell);
            elseif  size(self.workspace.zval, 1) == self.npnt && size(self.workspace.zval, 2) == 1
                % The corresponding points on all ellipsoids have the same z values.
                self.workspace.zval = repmat(self.workspace.zval, 1, self.nell);
            elseif  size(self.workspace.zval, 1) == 1 && size(self.workspace.zval, 2) == self.nell
                % All points on each ellipse have the same z values.
                self.workspace.zval = repmat(self.workspace.zval, self.npnt, 1);
            elseif  size(self.workspace.zval, 1) ~= self.npnt || size(self.workspace.zval, 2) ~= self.nell
                help("pm.vis.subplot.Ellipse3")
                disp("size(self.workspace.zval)")
                disp( size(self.workspace.zval) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``zval`` is incompatible with" + newline ...
                        + "the inferred number of ellipsoids ``nell`` in the input self.workspace." + newline ...
                        + "For more information, see the documentation of the ``zval`` component displayed above." + newline ...
                        + newline ...
                        );
            end

           %%%%
           %%%% Check the consistency of shape of ``self.ellindex``.
           %%%%

            if ~isempty(self.ellindex)
                self.workspace.ellindex = self.ellindex;
            else
                if ~pm.array.len(self.axes.colorScale)
                    self.axes.colorScale = "log";
                end
                self.workspace.ellindex = pm.array.logrange(1, self.nell, 150);
            end

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
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
            %       pm.vis.subplot.Ellipse3.make(varargin);
            %
            %   Example
            %   -------
            %
            %       p = pm.vis.subplot.Ellipse3();
            %       p.make("npnt", 200)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self.premake(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout`` and ``dfref``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%
            %%%% Generate column indices.
            %%%%

            fields = ["colx", "coly", "colz", "colc"];
            for icol = 1 : length(fields)
                field = fields(icol);
                self.(field) = length(fields) * ((1 : length(self.workspace.ellindex)) - 1) + icol;
            end
            %zmin = min(min(self.workspace.zval(:, self.colz)));
            %zmax = max(min(self.workspace.zval(:, self.colz)));
            %zrange = zmax - zmin;

            %%%%
            %%%% Generate ellipse self.workspace.
            %%%%

            %ellindexRange = log(self.workspace.ellindex(end)) - log(self.workspace.ellindex(1));
            bcrd = zeros(size(self.workspace.zval, 1), length(fields) * length(self.workspace.ellindex));
            for i = 1 : length(self.workspace.ellindex)
                iell = self.workspace.ellindex(i);
                icol = (i - 1) * length(fields) + 1;
                bcrd(:, icol : icol + 1) = pm.geom.ell2.getBorder(self.workspace.gramian(self.dims, self.dims, iell), self.workspace.center(self.dims, iell), size(self.workspace.zval, 1));
                bcrd(:, icol + 2) = self.workspace.zval(:, iell);
                if  3 < length(fields)
                    bcrd(:, icol + 3) = bcrd(:, icol + 2); %zmin + zrange * (log(iell)^3 - self.workspace.ellindex(1)^3) / ellindexRange;
                end
            end

            %%%%
            %%%% Generate the dataframe.
            %%%%

            self.df = pm.data.DataFrame(bcrd);

            %%%%
            %%%% Set the visualization axes names.
            %%%%

            if ~isempty(self.names)
                self.workspace.names = string(self.names);
                self.workspace.names = self.workspace.names(:);
            else
                self.workspace.names = strings(self.ndim, 1);
                for idim = 1 : self.ndim
                    self.workspace.names(idim) = "Dimension " + string(idim);
                end
            end

            if ~pm.array.len(self.xlabel.txt)
                self.xlabel.txt = self.workspace.names(self.dims(1));
            end
            if ~pm.array.len(self.ylabel.txt)
                self.ylabel.txt = self.workspace.names(self.dims(2));
            end
            if ~pm.array.len(self.zlabel.txt)
                self.zlabel.txt = "Z";
            end

            %%%%
            %%%% Set the visualization coloring.
            %%%%

            %if  isempty(self.colormap.enabled) || self.colormap.enabled
            %    if  pm.introspection.istype(self.colormap.cmap, "string")
            %        try
            %            cmap = eval(string(self.colormap.cmap)+"(length(self.ellindex))");
            %        catch
            %            error   ( "Failed to generate the color-mapping given the requested colormap: " + self.colormap.cmap ...
            %                    + "Please make sure the colormap string value is colormapping name recognized by the MATLAB's " ...
            %                    + "colormap() function." ...
            %                    );
            %        end
            %    else
            %        cmapsize = size(self.colormap.cmap);
            %        if isnumeric(self.colormap.cmap) && length(cmapsize)==2 && cmapsize(1)==rowindexLen
            %            cmap = self.colormap.cmap;
            %        else
            %            error   ( "A numeric value for the colormap must be given in the form of GRB triplet matrix, " ...
            %                    + "whose number of rows is the number of rows of the input dataframe to be visualized and, " ...
            %                    + "whose columns represent an RGB triplet." ...
            %                    );
            %        end
            %    end
            %end

            %%%%
            %%%% Generate the visualization.
            %%%%

            make@pm.vis.subplot.LineScatter3(self);
            if true % cval is empty
                self.fout.colorbar.Ticks = self.fout.axes.ZTick;
                %self.fout.colorbar.Limits = self.fout.axes.ZLim;
            end

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end