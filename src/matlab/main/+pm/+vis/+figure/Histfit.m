classdef Histfit < pm.vis.axes.BaseDF
    %
    %   This is the Histfit class for generating
    %   instances of 2-dimensional Histfit plots
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
    %       An object of ``pm.vis.axes.Histfit`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.Histfit(dfref);
    %       p = pm.vis.axes.Histfit(dfref, []);
    %       p = pm.vis.axes.Histfit(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Histfit(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.axes.BaseDF("Histfit", dfref);
        end
    end
end