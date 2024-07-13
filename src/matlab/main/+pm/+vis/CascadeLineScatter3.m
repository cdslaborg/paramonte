%>  \brief
%>  This is the CascadeLineScatter3 class for generating
%>  instances of 3-dimensional LineScatter3 [Cascade visualizations](@ref Cascade)
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
classdef CascadeLineScatter3 < pm.vis.Cascade

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.CascadeLineScatter3](@ref CascadeLineScatter3).<br>
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the class [pm.vis.Plot](@ref Plot).<br>
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.CascadeLineScatter3](@ref CascadeLineScatter3).<br>
        %>
        %>  \interface{CascadeLineScatter3}
        %>  \code{.m}
        %>
        %>      p = pm.vis.CascadeLineScatter3(dfref);
        %>      p = pm.vis.CascadeLineScatter3(dfref, varargin);
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
        %>  \example{CascadeLineScatter3}
        %>  \include{lineno} example/vis/CascadeLineScatter3/main.m
        %>  \vis{CascadeLineScatter3}
        %>  \image html example/vis/CascadeLineScatter3/CascadeLineScatter3.window.1.png width=700
        %>  \image html example/vis/CascadeLineScatter3/CascadeLineScatter3.window.2.png width=700
        %>  \image html example/vis/CascadeLineScatter3/CascadeLineScatter3.window.3.png width=700
        %>
        %>  \final{CascadeLineScatter3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:36 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = CascadeLineScatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Cascade(pm.vis.PlotLineScatter3(dfref), varargin{:});
        end

    end

end