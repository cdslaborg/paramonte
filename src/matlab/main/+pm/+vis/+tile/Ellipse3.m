classdef Ellipse3 < pm.vis.tile.Isotile
    %
    %   This is the Ellipse3 class for generating
    %   instances of 3-dimensional Ellipse3 tiles
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       gramian
    %
    %           See the corresponding input argument to the class ``pm.vis.subplot.Ellipse3``.
    %
    %       center
    %
    %           See the corresponding input argument to the class ``pm.vis.subplot.Ellipse3``.
    %
    %       zval
    %
    %           See the corresponding input argument to the class ``pm.vis.subplot.Ellipse3``.
    %
    %       cval
    %
    %           See the corresponding input argument to the class ``pm.vis.subplot.Ellipse3``.
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the parent object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``make()`` method.
    %
    %       \note
    %
    %           The input ``varargin`` can also contain the components
    %           of the ``subplot`` component of the parent object.
    %
    %   Attributes
    %   ----------
    %
    %       See below and also the documentation of the
    %       attributes of the superclass ``pm.vis.figure.Figure``.
    %
    %   Returns
    %   -------
    %
    %       An object of class ``pm.vis.tile.Ellipse3``.
    %
    %   Interface
    %   ---------
    %
    %       t = pm.vis.tile.Ellipse3();
    %       t = pm.vis.tile.Ellipse3(gramian);
    %       t = pm.vis.tile.Ellipse3(gramian, center);
    %       t = pm.vis.tile.Ellipse3(gramian, center, zval);
    %       t = pm.vis.tile.Ellipse3(gramian, center, zval, cval);
    %       t = pm.vis.tile.Ellipse3(gramian, center, zval, cval, varargin);
    %
    %   Example
    %   -------
    %
    %       t = pm.vis.tile.Ellipse3();
    %       t.make("dims", [1, 2]);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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
            self = self@pm.vis.tile.Isotile(pm.vis.subplot.Ellipse3(gramian, center, zval, cval), varargin{:});
        end
    end
end