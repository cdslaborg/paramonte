%>  \brief
%>  This is the SubplotContour3 class for generating
%>  instances of 3-dimensional Contour [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotContour3](@ref SubplotContour3::SubplotContour3) for example usage.<br>
%>
%>  \note
%>  See the documentation of the attributes
%>  of the superclass [pm.vis.Subplot](@ref Subplot).<br>
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
classdef SubplotContour3 < pm.vis.Subplot

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotContour3](@ref SubplotContour3).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.SubplotContour3](@ref SubplotContour3).<br>
        %>
        %>  \param[in]      dfref       :   See the documentation of the corresponding input
        %>                                  argument of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                    :   The output object of class [pm.vis.SubplotContour3](@ref SubplotContour3).<br>
        %>
        %>  \interface{SubplotContour3}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotContour3(dfref);
        %>      s = pm.vis.SubplotContour3(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \example{SubplotContour3}
        %>  \include{lineno} example/vis/SubplotContour3/main.m
        %>  \vis{SubplotContour3}
        %>  <br>\image html example/vis/SubplotContour3/SubplotContour3.1.png width=700
        %>
        %>  \final{SubplotContour3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotContour3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Subplot("Contour3", dfref, varargin{:});
        end

    end

end