classdef LineScatter < pm.vis.subplot.Subplot
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
    %       An object of ``pm.vis.subplot.LineScatter`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.LineScatter(dfref);
    %       p = pm.vis.subplot.LineScatter(dfref, varargin);
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
            self = self@pm.vis.subplot.Subplot("LineScatter", dfref, varargin{:});
        end
    end
end