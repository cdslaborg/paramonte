%>  \brief
%>  This is the PlotHeatmap class for generating
%>  instances of 2-dimensional Heatmap [Plot visualizations](@ref Plot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.PlotHeatmap](@ref PlotHeatmap::PlotHeatmap) for example usage.<br>
%>
%>  \note
%>  See the documentation of the attributes
%>  of the superclass [pm.vis.Plot](@ref Plot).<br>
%>
%>  \see
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Triplex](@ref Triplex)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef PlotHeatmap < pm.vis.Plot

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.PlotHeatmap](@ref PlotHeatmap).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.PlotHeatmap](@ref PlotHeatmap).<br>
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the class [pm.vis.Subplot](@ref Subplot).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output object of class [pm.vis.PlotHeatmap](@ref PlotHeatmap).<br>
        %>
        %>  \return
        %>  An object of [pm.vis.PlotHeatmap](@ref PlotHeatmap) class.<br>
        %>
        %>  \interface{PlotHeatmap}
        %>  \code{.m}
        %>
        %>      p = pm.vis.PlotHeatmap(dfref);
        %>      p = pm.vis.PlotHeatmap(dfref, varargin);
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
        %>  \example{PlotHeatmap}
        %>  \include{lineno} example/vis/PlotHeatmap/main.m
        %>  \vis{PlotHeatmap}
        %>  <br>\image html example/vis/PlotHeatmap/PlotHeatmap.1.png width=700
        %>
        %>  \final{PlotHeatmap}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:57 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = PlotHeatmap(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Plot(pm.vis.SubplotHeatmap(dfref), varargin{:});
        end

    end

end