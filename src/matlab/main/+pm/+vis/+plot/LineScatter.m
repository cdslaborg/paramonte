classdef LineScatter < pm.vis.plot.Plot
    %
    %   This is the LineScatter class for generating
    %   instances of 2-dimensional LineScatter plots
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
    %       An object of ``pm.vis.subplot.LineScatter`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.LineScatter(dfref);
    %       s = pm.vis.subplot.LineScatter(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = LineScatter(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.plot.Plot(pm.vis.subplot.LineScatter(dfref), varargin{:});
        end
    end
end