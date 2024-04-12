classdef Histogram2 < pm.vis.AxesData
    %
    %   This is the Histogram2 class for generating
    %   instances of 3-dimensional Histogram2 plots
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
    %       An object of ``pm.vis.axes.Histogram2`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Histogram2(dfref);
    %       p = pm.vis.axes.Histogram2(dfref, []);
    %       p = pm.vis.axes.Histogram2(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Histogram2(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.AxesData("Histogram2", dfref);
        end
    end
end