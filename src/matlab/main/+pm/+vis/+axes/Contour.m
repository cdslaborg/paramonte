classdef Contour < pm.vis.axes.BaseDF
    %
    %   This is the Contour class for generating
    %   instances of 2-dimensional Contour plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       dfref
    %
    %           See the documentation of the corresponding input
    %           argument of the parent class ``pm.vis.axes.BaseDF``.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the parent class ``pm.vis.axes.BaseDF``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.axes.Contour`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Contour(dfref);
    %       p = pm.vis.axes.Contour(dfref, []);
    %       p = pm.vis.axes.Contour(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Contour(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.axes.BaseDF("Contour", dfref);
        end
    end
end