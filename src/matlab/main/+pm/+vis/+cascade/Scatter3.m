%
%   This is the Scatter3 class for generating
%   instances of 3-dimensional Scatter3 plots
%   based on the relevant MATLAB
%   intrinsic functions.
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
%>  \return
%       An object of ``pm.vis.cascade.Scatter3`` class.
%
%   Interface
%   ---------
%
%       s = pm.vis.cascade.Scatter3(dfref);
%       s = pm.vis.cascade.Scatter3(dfref, varargin);
%
%>  \final{}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Scatter3 < pm.vis.cascade.Cascade
    methods(Access = public)
        function self = Scatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.cascade.Cascade(pm.vis.plot.Scatter3(dfref), varargin{:});
        end
    end
end