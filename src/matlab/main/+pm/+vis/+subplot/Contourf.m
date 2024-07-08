%>  \brief
%>  This is the Contourf class for generating
%>  instances of 2-dimensional Contourf [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.subplot.Contourf](@ref Contourf::Contourf) for example usage.<br>
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
classdef Contourf < pm.vis.subplot.Subplot
    methods(Access = public)
        %>  \brief
        %>  Construct and return an object of class [pm.vis.subplot.Contourf](@ref Contourf).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.subplot.Contourf](@ref Contourf).<br>
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
        %>  An object of [pm.vis.subplot.Contourf](@ref Contourf) class.<br>
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.subplot.Subplot](@ref Subplot).<br>
        %>
        %>  \interface{Contourf}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Contourf(dfref);
        %>      s = pm.vis.subplot.Contourf(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Contourf}
        %>  \include{lineno} example/vis/subplot/Contourf/main.m
        %>  \vis{Subplot}
        %>  <br>\image html example/vis/subplot/Contourf/Subplot.Contourf.1.png width=700
        %>
        %>  \final{Contourf}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Contourf(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Contourf", dfref, varargin{:});
        end
    end
end