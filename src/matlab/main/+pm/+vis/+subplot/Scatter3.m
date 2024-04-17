classdef Scatter3 < pm.vis.subplot.Subplot
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
    %           argument of the parent class ``pm.vis.subplot.Subplot``.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the parent class ``pm.vis.subplot.Subplot``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.subplot.Scatter3`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.Scatter3(dfref);
    %       p = pm.vis.subplot.Scatter3(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Scatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Scatter3", dfref, varargin{:});
        end
    end
end