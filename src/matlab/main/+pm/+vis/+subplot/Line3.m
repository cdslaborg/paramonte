classdef Line3 < pm.vis.subplot.Subplot
    %
    %   This is the Line3 class for generating
    %   instances of 3-dimensional Line3 plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       dfref
    %
    %           See the documentation of the corresponding input
    %           argument of the superclass ``pm.vis.subplot.Subplot``.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the superclass ``pm.vis.subplot.Subplot``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.subplot.Line3`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Line3(dfref);
    %       s = pm.vis.subplot.Line3(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Line3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Line3", dfref, varargin{:});
        end
    end
end