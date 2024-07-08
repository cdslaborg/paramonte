%>  \brief
%>  This is the LineScatter3 class for generating
%>  instances of 3-dimensional Line-Scatter [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.subplot.LineScatter3](@ref LineScatter3::LineScatter3) for example usage.<br>
%>
%>  \see
%>  [pm.vis.cascade](@ref \psldir/main/+pm/+vis/+cascade)<br>
%>  [pm.vis.subplot](@ref \psldir/main/+pm/+vis/+subplot)<br>
%>  [pm.vis.figure](@ref \psldir/main/+pm/+vis/+figure)<br>
%>  [pm.vis.corner](@ref \psldir/main/+pm/+vis/+corner)<br>
%>  [pm.vis.plot](@ref \psldir/main/+pm/+vis/+plot)<br>
%>  [pm.vis.tile](@ref \psldir/main/+pm/+vis/+tile)<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef LineScatter3 < pm.vis.subplot.Subplot
    methods(Access = public)
        %>  \brief
        %>  Construct and return an object of class [pm.vis.subplot.LineScatter3](@ref LineScatter3).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.subplot.LineScatter3](@ref LineScatter3).<br>
        %>
        %>  \param[in]      dfref       :   See the documentation of the corresponding input
        %>                                  argument of the superclass [pm.vis.subplot.Subplot](@ref Subplot).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  An object of [pm.vis.subplot.LineScatter3](@ref LineScatter3) class.<br>
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.subplot.Subplot](@ref Subplot).<br>
        %>
        %>  \interface{LineScatter3}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.LineScatter3(dfref);
        %>      s = pm.vis.subplot.LineScatter3(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{LineScatter3}
        %>  \include{lineno} example/vis/subplot/LineScatter3/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/subplot/LineScatter3/Subplot.LineScatter3.1.png width=700
        %>
        %>  \final{LineScatter3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = LineScatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("LineScatter3", dfref, varargin{:});
        end
    end
end