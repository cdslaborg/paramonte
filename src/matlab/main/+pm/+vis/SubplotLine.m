%>  \brief
%>  This is the SubplotLine class for generating
%>  instances of 2-dimensional Line [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotLine](@ref SubplotLine::SubplotLine) for example usage.<br>
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
classdef SubplotLine < pm.vis.Subplot

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotLine](@ref SubplotLine).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.SubplotLine](@ref SubplotLine).<br>
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
        %>  ``self``                    :   The output object of class [pm.vis.SubplotLine](@ref SubplotLine).<br>
        %>
        %>  \interface{SubplotLine}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotLine(dfref);
        %>      s = pm.vis.SubplotLine(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \example{SubplotLine}
        %>  \include{lineno} example/vis/SubplotLine/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/SubplotLine/SubplotLine.1.png width=700
        %>
        %>  \final{SubplotLine}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotLine(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Subplot("Line", dfref, varargin{:});
        end

    end

end