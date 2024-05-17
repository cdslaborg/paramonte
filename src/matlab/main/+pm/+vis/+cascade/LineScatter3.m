classdef LineScatter3 < pm.vis.cascade.Cascade
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
    %           argument of the class ``pm.vis.plot.Plot``.
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the parent object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``make()`` method.
    %
    %       \note
    %
    %           The input ``varargin`` can also contain the components
    %           of the ``template`` component of the parent object.
    %
    %   Attributes
    %   ----------
    %
    %       See the documentation of the attributes
    %       of the superclass ``pm.vis.cascade.Cascade``.
    %
    %   Returns
    %   -------
    %
    %       An object of ``pm.vis.cascade.LineScatter3`` class.
    %
    %   Interface
    %   ---------
    %
    %       p = pm.vis.cascade.LineScatter3(dfref);
    %       p = pm.vis.cascade.LineScatter3(dfref, varargin);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = LineScatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.cascade.Cascade(pm.vis.plot.LineScatter3(dfref), varargin{:});
        end
    end
end