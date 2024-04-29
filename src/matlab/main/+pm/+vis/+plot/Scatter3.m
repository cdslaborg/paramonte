classdef Scatter3 < pm.vis.plot.Plot
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
    %           argument of the superclass ``pm.vis.plot.Plot``.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the superclass ``pm.vis.plot.Plot``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.subplot.Scatter3`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Scatter3(dfref);
    %       s = pm.vis.subplot.Scatter3(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Scatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.plot.Plot(pm.vis.subplot.Scatter3(dfref), varargin{:});
        end
    end
end