classdef Contourf < pm.vis.subplot.Subplot
    %
    %   This is the Contourf class for generating
    %   instances of 2-dimensional Contourf plots
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
    %       An object of ``pm.vis.subplot.Contourf`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.subplot.Contourf(dfref);
    %       p = pm.vis.subplot.Contourf(dfref, []);
    %       p = pm.vis.subplot.Contourf(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = Contourf(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Contourf", dfref);
        end
    end
end