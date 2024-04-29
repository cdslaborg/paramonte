classdef Histogram < pm.vis.subplot.Subplot
    %
    %   This is the Histogram class for generating
    %   instances of 2-dimensional Histogram plots
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
    %       An object of ``pm.vis.subplot.Histogram`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Histogram(dfref);
    %       s = pm.vis.subplot.Histogram(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Histogram(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Histogram", dfref, varargin{:});
        end
    end
end