classdef Contourf < pm.vis.AxesData
    %
    %   This is the Contourf class for generating
    %   instances of 2-dimensional Contourf plots
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
    %       An object of ``pm.vis.axes.Contourf`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Contourf(dfref);
    %       p = pm.vis.axes.Contourf(dfref, []);
    %       p = pm.vis.axes.Contourf(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Contourf(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.AxesData("Contourf", dfref);
        end
    end
end