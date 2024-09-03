%>  \brief
%>  This is the CascadeContourf class for generating
%>  instances of 2-dimensional Contourf [Cascade visualizations](@ref Cascade)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.vis.Cascade](@ref Cascade).<br>
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
classdef CascadeContourf < pm.vis.Cascade

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.CascadeContourf](@ref CascadeContourf).<br>
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
        %>  ``self``                :   The output scalar object of class [pm.vis.CascadeContourf](@ref CascadeContourf).<br>
        %>
        %>  \interface{CascadeContourf}
        %>  \code{.m}
        %>
        %>      p = pm.vis.CascadeContourf(dfref);
        %>      p = pm.vis.CascadeContourf(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``template`` component of the parent object.<br>
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.Cascade](@ref Cascade).<br>
        %>
        %>  \example{CascadeContourf}
        %>  \include{lineno} example/vis/CascadeContourf/main.m
        %>  \vis{CascadeContourf}
        %>  \image html example/vis/CascadeContourf/CascadeContourf.window.1.png width=700
        %>  \image html example/vis/CascadeContourf/CascadeContourf.window.2.png width=700
        %>  \image html example/vis/CascadeContourf/CascadeContourf.window.3.png width=700
        %>
        %>  \final{CascadeContourf}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:16 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = CascadeContourf(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Cascade(pm.vis.PlotContourf(dfref), varargin{:});
        end

    end

end