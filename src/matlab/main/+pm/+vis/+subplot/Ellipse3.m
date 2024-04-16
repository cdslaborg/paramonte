classdef Ellipse3 < pm.vis.axes.LineScatter3
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
    %           The MATLAB (function handle returning a)
    %           3D rectangle of shape ``[ndim, ndim, nell]``, each
    %           subset ``gramian(1:ndim, 1:ndim, igram)`` of which represents
    %           the ``igram`` Gramian matrix of an 2D-planar ellipse to visualize,
    %           where ``ndim`` refers to the number of dimensions of the
    %           hyper-ellipsoid represented by the ``igram`` Gramian
    %           and ``nell`` is the number of Gramian matrices.
    %
    %           \note
    %
    %               While it is possible to pass a 3D rectangle directly to this class,
    %               it is highly recommended to pass a function handle that returns
    %               such data when called, allowing the visualization data to be
    %               dynamically updated when needed.
    %
    %       center
    %
    %           The input MATLAB (function handle returning a)
    %           vector of doubles of size ``ndim`` or matrix of
    %           shape ``[ndim, nell]`` containing the 2D coordinates
    %           of the centers of the 2D ellipsoids to display,
    %           where ``ndim`` refers to the number of dimensions of the
    %           hyper-ellipsoid represented by the ``igram`` Gramian
    %           and ``nell`` is the number of Gramian matrices.
    %
    %           \note
    %
    %               While it is possible to pass a 2D rectangle directly to this class,
    %               it is highly recommended to pass a function handle that returns
    %               such data when called, allowing the visualization data to be
    %               dynamically updated when needed.
    %
    %       zval
    %
    %           The input MATLAB (function handle returning a)
    %           vector doubles of size ``nell`` or matrix of
    %           shape ``[nell]`` containing the z-axis coordinates
    %           of each of the 2D ellipsoids to display, where ``nell``
    %           is the number of Gramian matrices in the input ``gramian``.
    %           (**optional**. The default is the values within the ``indices`` component of the object.)
    %
    %           \note
    %
    %               While it is possible to pass a 1D vector directly to this class,
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
    %       of the parent class ``pm.vis.axes.LineScatter3``.
    %
    %   Returns
    %   -------
    %
    %       An object of class ``pm.vis.axes.Ellipse3``.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Ellipse3(gramian, center);
    %       p = pm.vis.axes.Ellipse3(gramian, center, varargin);
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
        %
        npnt = 50;
        %
        %   dims
        %
        %       The vector of size ``2`` of MATLAB whole numbers representing the dimensions
        %       of the input ``gramian`` and ``center`` data to use for 2D ellipsoid visualization.
        %       The default is ``[1, 2]``.
        %
        dims = [1, 2];
        %
        %   indices
        %
        %       The vector of MATLAB whole numbers representing the indices of the 2D
        %       ellipsoids to display, represented by the input ``gramian`` and ``center``.
        %       The specified indices will serve as the visualization values on the z-axis.
        %       The default is ``indices = pm.array.logrange(start = 1, stop = nell, count = 100);``,
        %       where ``nell`` represents the number of ellipsoids to visualize.
        %
        indices = [];
        %
        %   names
        %
        %       The vector of MATLAB strings containing the names of dimensions
        %       of the space within which the Gramian matrices are defined.
        %       The default is ``"Dimension i"`` where ``i`` is replaced
        %       by the ID of the corresponding axis.
        %
        names = [];
    end

    properties (Access = Hidden)
        gramian = [];
        center = [];
        zval = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetint(self, varargin)
            self.indices = struct();
            self.dims = [1, 2];
            self.npnt = 50;
            self.names = [];
            resetint@pm.vis.LineScatter3(self, varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Ellipse3(gramian, center, zval, varargin)
            varargin = {varargin{:}, "gramian", gramian, "center", center, "zval", zval};
            self = self@pm.vis.axes.LineScatter3([], varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
            %
            %   Generate a plot from the selected columns
            %   of the parent object component ``dfref``.
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
            if self.isdryrun;
                make@pm.vis.LineScatter3(self, varargin{:});
                return;
            elseif ~isempty(varargin)
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%% Set the data.

            if  isa(self.gramian, 'function_handle')
                gramian = self.gramian();
            else
                gramian = self.gramian;
            end
            if  isa(self.center, 'function_handle')
                center = self.center();
            else
                center = self.center;
            end
            if  isa(self.zval, 'function_handle')
                zval = self.zval();
            else
                zval = self.zval;
            end

            %%%% Ensure the correctness of the user-specified data.

            asserted = ~isempty(gramian) && ~isempty(center);
            if  asserted
                asserted = asserted && size(gramian, 1) == size(center, 1);
                asserted = asserted && size(gramian, 2) == size(center, 1);
                asserted = asserted &&(size(gramian, 3) == size(center, 2) || size(gramian, 3) == 1 || size(center, 2) == 1);
            end

            if ~asserted
                help("pm.vis.axes.Ellipse3")
                disp("size(gramian)")
                disp( size(gramian) )
                disp("size(center)")
                disp( size(center) )
                error   ( newline ...
                        + "The shapes of the specified ``gramian`` and ``center`` are incompatible." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            else
                ndim = size(gramian, 2);
                nell = max(size(gramian, 3), size(center, 2));
            end

           %%%% Check the minimum size of ``size(gramian, 1)`` and ``size(center, 1)``.

            if  length(self.dims) ~= 2
                help("pm.vis.axes.Ellipse3")
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

            %%%% Set up dimension names.

            if ~isempty(self.names)
                names = self.names;
                if  length(names) ~= ndim
                    help("pm.vis.axes.Ellipse3")
                    disp("self.names")
                    disp( self.names )
                    disp("size(self.names)")
                    disp( size(self.names) )
                    error   ( newline ...
                            + "The component ``names`` must be a vector of size ``ndim`` of strings." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end
            else
                names = strings(ndim, 1);
                for idim = 1, ndim
                    self.names(idim) = "Dimension " + string(idim);
                end
            end

            %%%% Set up dimension names.

            make@pm.vis.LineScatter3(self);

        end

    end

end