%>  \brief
%>  This is the Ellipse class for generating
%>  instances of 2-dimensional Ellipse plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Ellipse < pm.vis.plot.Plot
    methods(Access = public)
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class [pm.vis.subplot.Ellipse](@ref Ellipse).
        %>
        %>  \param[in]  center      :   See the corresponding input argument to the class [pm.vis.subplot.Ellipse](@ref Ellipse).
        %>
        %>  \param[in]  cval        :   See the corresponding input argument to the class [pm.vis.subplot.Ellipse](@ref Ellipse).
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \note
        %>  See below and also the documentation of the
        %>  attributes of the superclass [pm.vis.figure.Figure](@ref Figure).
        %>
        %>  \return
        %>  An object of class [pm.vis.plot.Ellipse](@ref Ellipse).
        %>
        %>  \interface{Ellipse}
        %>  \code{.m}
        %>
        %>      p = pm.vis.plot.Ellipse();
        %>      p = pm.vis.plot.Ellipse(gramian);
        %>      p = pm.vis.plot.Ellipse(gramian, center);
        %>      p = pm.vis.plot.Ellipse(gramian, center, cval);
        %>      p = pm.vis.plot.Ellipse(gramian, center, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Ellipse}
        %>
        %>      p = pm.vis.plot.Ellipse();
        %>      p.make("dims", [1, 2]);
        %>
        %>  \final{Ellipse}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:53 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Ellipse(gramian, center, cval, varargin)
            %%%% Define the missing optional values as empty with the right rank.
            if  nargin < 3
                cval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end
            self = self@pm.vis.plot.Plot(pm.vis.subplot.Ellipse(gramian, center, cval), varargin{:});
        end
    end
end