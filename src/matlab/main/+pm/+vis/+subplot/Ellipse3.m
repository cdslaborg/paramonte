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
    %           2D rectangle of shape ``[ndim, ndim]`` or
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
    %           vector of MATLAB doubles of size ``ndim`` or matrix of
    %           shape ``[ndim, nell]`` containing the 2D coordinates of the
    %           centers of the 2D ellipsoids to display.
    %
    %           \note
    %
    %               While it is possible to pass a 2D rectangle directly to this class,
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
        %   indices
        %
        %       The vector of MATLAB whole numbers representing the indices of the 2D
        %       ellipsoids to display, represented by the input ``gramian`` and ``center``.
        %       The default is ``indices = pm.array.logrange(start = 1, stop = nell, count = 100);``.
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
        gramref = [];
        centref = [];
    end

    methods (Access = public)

        function self = Ellipse3(gramian, center, varargin)

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
            self = self@pm.vis.axes.LineScatter3([], varargin{:});
        end

    end

end