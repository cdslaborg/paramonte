%>  \brief
%>  This is the SubplotContourf class for generating
%>  instances of 2-dimensional Contour [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotContourf](@ref SubplotContourf::SubplotContourf) for example usage.<br>
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
classdef SubplotContourf < pm.vis.Subplot

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotContourf](@ref SubplotContourf).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.SubplotContourf](@ref SubplotContourf).<br>
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
        %>  ``self``                    :   The output object of class [pm.vis.SubplotContourf](@ref SubplotContourf).<br>
        %>
        %>  \interface{SubplotContourf}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotContourf(dfref);
        %>      s = pm.vis.SubplotContourf(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \example{SubplotContourf}
        %>  \include{lineno} example/vis/SubplotContourf/main.m
        %>  \vis{SubplotContourf}
        %>  <br>\image html example/vis/SubplotContourf/SubplotContourf.1.png width=700
        %>
        %>  \final{SubplotContourf}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotContourf(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Subplot("Contourf", dfref, varargin{:});
        end

    end

end