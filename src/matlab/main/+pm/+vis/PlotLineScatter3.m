%>  \brief
%>  This is the PlotLineScatter3 class for generating
%>  instances of 3-dimensional Line-Scatter [Plot visualizations](@ref Plot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.PlotLineScatter3](@ref PlotLineScatter3::PlotLineScatter3) for example usage.<br>
%>
%>  \note
%>  See the documentation of the attributes
%>  of the superclass [pm.vis.Plot](@ref Plot).<br>
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
%>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef PlotLineScatter3 < pm.vis.Plot

    methods(Access = public)

        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the superclass [pm.vis.Plot](@ref Plot).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output object of class [pm.vis.PlotLineScatter3](@ref PlotLineScatter3).<br>
        %>
        %>  \interface{PlotLineScatter3}
        %>  \code{.m}
        %>
        %>      p = pm.vis.PlotLineScatter3(dfref);
        %>      p = pm.vis.PlotLineScatter3(dfref, varargin);
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
        %>  \example{PlotLineScatter3}
        %>  \include{lineno} example/vis/PlotLineScatter3/main.m
        %>  \vis{PlotLineScatter3}
        %>  <br>\image html example/vis/PlotLineScatter3/PlotLineScatter3.1.png width=700
        %>
        %>  \final{PlotLineScatter3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:12 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = PlotLineScatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Plot(pm.vis.SubplotLineScatter3(dfref), varargin{:});
        end

    end

end