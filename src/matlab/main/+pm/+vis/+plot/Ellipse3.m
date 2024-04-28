classdef Ellipse3 < pm.vis.subplot.LineScatter
    %
    %   This is the Ellipse3 class for generating
    %   instances of 3-dimensional Ellipse3 plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       gramref
    %
    %           The input MATLAB array that can be:
    %
    %               A (function handle returning a) 3D rectangle of shape ``[ndim, ndim, ngram]``,
    %               each subset ``gramref(1:ndim, 1:ndim, igram)`` of which represents
    %               the ``igram`` Gramian matrix of an 2D-planar ellipse to visualize,
    %               where ``ndim`` refers to the number of dimensions of the
    %               hyper-ellipsoid represented by the ``igram`` Gramian
    %               and ``ngram`` is the number of Gramian matrices.
    %
    %           While it is possible to pass the 3D rectangle directly to this class,
    %           it is highly recommended to pass a function handle that returns
    %           such data when called, allowing the visualization data to be
    %           dynamically updated when needed.
    %
    %   Attributes
    %   ----------
    %
    %       See below and the documentation of the attributes
    %       of the parent class ``pm.vis.subplot.LineScatter``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``Ellipse3`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.Ellipse3(gramref);
    %       p = pm.vis.subplot.Ellipse3(gramref, []);
    %       p = pm.vis.subplot.Ellipse3(gramref, [], []);
    %       p = pm.vis.subplot.Ellipse3(gramref, [], zdfref);
    %       p = pm.vis.subplot.Ellipse3(gramref, meanref, zdfref);
    %       p = pm.vis.subplot.Ellipse3([], meanref, zdfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties(Access = public)
        %
        %   resolution
        %
        %       The optional scalar MATLAB whole-number representing
        %       the number of points used for visualization of ellipses.
        %       Higher resolution is particularly helpful with visualizing
        %       elongated ellipses. The default value is ``100``.
        %
        resolution = [];
    end

    properties(Access = Hidden)
        gramref = [];
        meanref = [];
        zdfref = [];
    end

    methods(Access = public)
        function self = Ellipse3(gramref, meanref, zdfref)
            if  nargin < 3
                zdfref = [];
            end
            if  nargin < 2
                meanref = [];
            end
            if  nargin < 1
                gramref = [];
            end
            if ~isempty(gramref)
                if ~isempty(gramref)
                if  size(gramref, 1) ~= size(gramref, 2)
                    help("pm.vis.subplot.Ellipse3")
                    error   ( newline ...
                            + "The specified ``gramref`` must be a 3D rectangle." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end
                if ~isempty(meanref)

                end
            && isempty(gramref) && isempty(zdfref)
            self = self@pm.vis.subplot.LineScatter(, gramref);
        end
    end
end