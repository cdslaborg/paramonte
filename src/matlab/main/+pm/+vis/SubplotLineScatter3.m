%>  \brief
%>  This is the SubplotLineScatter3 class for generating
%>  instances of 3-dimensional Line-Scatter [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3::SubplotLineScatter3) for example usage.<br>
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
classdef SubplotLineScatter3 < pm.vis.Subplot

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3).<br>
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
        %>  ``self``                    :   The output object of class [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3).<br>
        %>
        %>  \interface{SubplotLineScatter3}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotLineScatter3(dfref);
        %>      s = pm.vis.SubplotLineScatter3(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>
        %>  \example{SubplotLineScatter3}
        %>  \include{lineno} example/vis/SubplotLineScatter3/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/SubplotLineScatter3/SubplotLineScatter3.1.png width=700
        %>
        %>  \final{SubplotLineScatter3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotLineScatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Subplot("LineScatter3", dfref, varargin{:});
        end

    end

end