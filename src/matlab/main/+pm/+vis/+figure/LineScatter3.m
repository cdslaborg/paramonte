classdef LineScatter3 < pm.vis.axes.BaseDF
    %
    %   This is the LineScatter3 class for generating
    %   instances of 3-dimensional LineScatter3 plots
    %   based on the relevant MATLAB
    %   intrinsic functions.
    %
    %   Parameters
    %   ----------
    %
    %       dfref
    %
    %           See the documentation of the corresponding input
    %           argument of the parent class ``pm.vis.axes.BaseDF``.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the parent class ``pm.vis.axes.BaseDF``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.axes.LineScatter3`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.axes.LineScatter3(dfref);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods (Access = public)
        function self = LineScatter3(dfref)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.axes.BaseDF("LineScatter3", dfref);
        end
    end
end