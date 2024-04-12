classdef Scatter3 < pm.vis.AxesData
    %
    %   This is the Scatter3 class for generating
    %   instances of 3-dimensional Scatter3 plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       dfref
    %
    %           See the documentation of the corresponding input
    %           argument of the parent class ``pm.vis.AxesData``.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the parent class ``pm.vis.AxesData``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.axes.Scatter3`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Scatter3(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Scatter3(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.AxesData("Scatter3", dfref);
        end
    end
end