%>  \brief
%>  This is the PlotEllipse3 class for generating
%>  instances of 3-dimensional Ellipse [Plot visualizations](@ref Plot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.PlotEllipse3](@ref PlotEllipse3::PlotEllipse3) for example usage.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.vis.Plot](@ref Plot).<br>
%>
%>  \see
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Corner](@ref Corner)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef PlotEllipse3 < pm.vis.Plot

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.PlotEllipse3](@ref PlotEllipse3).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.PlotEllipse3](@ref PlotEllipse3).<br>
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>  \param[in]  center      :   See the corresponding input argument to the class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>  \param[in]  zval        :   See the corresponding input argument to the class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>  \param[in]  cval        :   See the corresponding input argument to the class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output object of class [pm.vis.PlotEllipse3](@ref PlotEllipse3).<br>
        %>
        %>  \interface{PlotEllipse3}
        %>  \code{.m}
        %>
        %>      p = pm.vis.PlotEllipse3();
        %>      p = pm.vis.PlotEllipse3(gramian);
        %>      p = pm.vis.PlotEllipse3(gramian, center);
        %>      p = pm.vis.PlotEllipse3(gramian, center, zval);
        %>      p = pm.vis.PlotEllipse3(gramian, center, zval, cval);
        %>      p = pm.vis.PlotEllipse3(gramian, center, zval, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.<br>
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.Plot](@ref Plot).<br>
        %>
        %>  \example{PlotEllipse3}
        %>  \include{lineno} example/vis/PlotEllipse3/main.m
        %>  \vis{PlotEllipse3}
        %>  <br>\image html example/vis/PlotEllipse3/PlotEllipse3.1.png width=700
        %>
        %>  \final{PlotEllipse3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:55 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = PlotEllipse3(gramian, center, zval, cval, varargin)
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
            self = self@pm.vis.Plot(pm.vis.SubplotEllipse3(gramian, center, zval, cval), varargin{:});
        end

    end

end