classdef Line < pm.vis.plot.Plot
    %
    %   This is the Line class for generating
    %   instances of 2-dimensional Line plots
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
    %       An object of ``pm.vis.subplot.Line`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Line(dfref);
    %       s = pm.vis.subplot.Line(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Line(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.plot.Plot(pm.vis.subplot.Line(dfref), varargin{:});
        end
    end
end