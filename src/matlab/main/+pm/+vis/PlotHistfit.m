%>  \brief
%>  This is the Histfit class for generating
%>  instances of 2-dimensional Histfit [Plot visualizations](@ref Plot)
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef PlotHistfit < pm.vis.Plot

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
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.Plot](@ref Plot).
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \return
        %>  An object of [pm.vis.PlotHistfit](@ref PlotHistfit) class.
        %>
        %>  \interface{PlotHistfit}
        %>  \code{.m}
        %>
        %>      p = pm.vis.PlotHistfit(dfref);
        %>      p = pm.vis.PlotHistfit(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{PlotHistfit}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:59 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = PlotHistfit(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Plot(pm.vis.SubplotHistfit(dfref), varargin{:});
        end

    end

end