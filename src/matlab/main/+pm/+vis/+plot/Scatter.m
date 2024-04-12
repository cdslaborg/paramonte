classdef Scatter < pm.vis.AxesData
    %
    %   This is the Scatter class for generating
    %   instances of 2-dimensional Scatter plots
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
    %       An object of ``Scatter`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Scatter(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Scatter(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.AxesData("Scatter", dfref);
        end
    end
end