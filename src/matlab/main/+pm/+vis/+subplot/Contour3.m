classdef Contour3 < pm.vis.subplot.Subplot
    %
    %   This is the Contour3 class for generating
    %   instances of 3-dimensional Contour3 plots
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
    %       An object of ``pm.vis.subplot.Contour3`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.Contour3(dfref);
    %       p = pm.vis.subplot.Contour3(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Contour3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Contour3", dfref, varargin{:});
        end
    end
end