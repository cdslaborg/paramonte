%>  \brief
%>  This is the CascadeHistogram2 class for generating
%>  instances of 2-dimensional Histogram2 [Cascade visualizations](@ref Cascade)
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
%>  [pm.vis.Triplex](@ref Triplex)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef CascadeHistogram2 < pm.vis.Cascade

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.CascadeHistogram2](@ref CascadeHistogram2).<br>
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the class [pm.vis.Subplot](@ref Subplot).<br>
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.CascadeHistogram2](@ref CascadeHistogram2).<br>
        %>
        %>  \interface{CascadeHistogram2}
        %>  \code{.m}
        %>
        %>      p = pm.vis.CascadeHistogram2(dfref);
        %>      p = pm.vis.CascadeHistogram2(dfref, varargin);
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
        %>  \example{CascadeHistogram2}
        %>  \include{lineno} example/vis/CascadeHistogram2/main.m
        %>  \vis{CascadeHistogram2}
        %>  \image html example/vis/CascadeHistogram2/CascadeHistogram2.window.1.png width=700
        %>  \image html example/vis/CascadeHistogram2/CascadeHistogram2.window.2.png width=700
        %>  \image html example/vis/CascadeHistogram2/CascadeHistogram2.window.3.png width=700
        %>
        %>  \final{CascadeHistogram2}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:31 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = CascadeHistogram2(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Cascade(pm.vis.PlotHistogram2(dfref), varargin{:});
        end

    end

end