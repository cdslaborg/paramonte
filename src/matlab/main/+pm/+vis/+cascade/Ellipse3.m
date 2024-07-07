%>  \brief
%>  This is the Ellipse3 class for generating
%>  instances of 3-dimensional Ellipse3 plots
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.vis.cascade.Cascade](@ref Cascade).<br>
%>
%>  \see
%>  [pm.vis.cascade](@ref \psldir/main/+pm/+vis/+cascade)<br>
%>  [pm.vis.subplot](@ref \psldir/main/+pm/+vis/+subplot)<br>
%>  [pm.vis.figure](@ref \psldir/main/+pm/+vis/+figure)<br>
%>  [pm.vis.corner](@ref \psldir/main/+pm/+vis/+corner)<br>
%>  [pm.vis.plot](@ref \psldir/main/+pm/+vis/+plot)<br>
%>  [pm.vis.tile](@ref \psldir/main/+pm/+vis/+tile)<br>
%>
%>  \final
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Ellipse3 < pm.vis.cascade.Cascade
    methods(Access = public)
        %>  \brief
        %>  Construct and return an object of class [pm.vis.cascade.Ellipse3](@ref Ellipse3).<br>
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class [pm.vis.plot.Ellipse3](@ref Ellipse3).<br>
        %>  \param[in]  center      :   See the corresponding input argument to the class [pm.vis.plot.Ellipse3](@ref Ellipse3).<br>
        %>  \param[in]  zval        :   See the corresponding input argument to the class [pm.vis.plot.Ellipse3](@ref Ellipse3).<br>
        %>  \param[in]  cval        :   See the corresponding input argument to the class [pm.vis.plot.Ellipse3](@ref Ellipse3).<br>
        %>  \param[in]  varargin    :    Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.cascade.Ellipse3](@ref Ellipse3).<br>
        %>
        %>  \interface{Ellipse3}
        %>  \code{.m}
        %>
        %>      p = pm.vis.cascade.Ellipse3();
        %>      p = pm.vis.cascade.Ellipse3(gramian);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center, zval);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center, zval, cval);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center, zval, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``template`` component of the parent object.<br>
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.cascade.Cascade](@ref Cascade).<br>
        %>
        %>  \example{Ellipse3}
        %>  \include{lineno} example/vis/cascade/Ellipse3/main.m
        %>  \vis{Ellipse3}
        %>  \image html example/vis/cascade/Ellipse3/Ellipse3.window.1.png width=700
        %>  \image html example/vis/cascade/Ellipse3/Ellipse3.window.2.png width=700
        %>  \image html example/vis/cascade/Ellipse3/Ellipse3.window.3.png width=700
        %>  \image html example/vis/cascade/Ellipse3/Ellipse3.window.4.png width=700
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:21 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Ellipse3(gramian, center, zval, cval, varargin)
            %%%% Define the missing optional values as empty with the right rank.
            if  nargin < 4
                cval = zeros(0, 0);
            end
            if  nargin < 3
                zval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end
            self = self@pm.vis.cascade.Cascade(pm.vis.plot.Ellipse3(gramian, center, zval, cval), varargin{:});
        end
    end
end