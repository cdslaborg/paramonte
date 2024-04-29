classdef Histfit < pm.vis.subplot.Subplot
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
    %       An object of ``pm.vis.subplot.Histfit`` class.
    %
    %   Interface
    %   ---------
    %
    %       s = pm.vis.subplot.Histfit(dfref);
    %       s = pm.vis.subplot.Histfit(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Histfit(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Histfit", dfref, varargin{:});
        end
    end
end