classdef Histogram2 < pm.vis.plot.Plot
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
    %       An object of ``pm.vis.subplot.Histogram2`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Histogram2(dfref);
    %       s = pm.vis.subplot.Histogram2(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Histogram2(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.plot.Plot(pm.vis.subplot.Histogram2(dfref), varargin{:});
        end
    end
end