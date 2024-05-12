classdef Contour3 < pm.vis.tile.Isotile
    %
    %   This is the Contour3 class for generating
    %   instances of 3-dimensional Contour3 tiles
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       dfref
    %
    %           See the documentation of the corresponding input
    %           argument of the superclass ``pm.vis.tile.Isotile``.
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
    %       See the documentation of the attributes
    %       of the superclass ``pm.vis.tile.Isotile``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.tile.Contour3`` class.
    %
    %   Interface
    %   ---------
    %
    %       t = pm.vis.tile.Contour3(dfref);
    %       t = pm.vis.tile.Contour3(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Contour3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.tile.Isotile(pm.vis.subplot.Contour3(dfref), varargin{:});
        end
    end
end