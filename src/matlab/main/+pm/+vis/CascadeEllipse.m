%>  \brief
%>  This is the CascadeEllipse class for generating
%>  instances of 2-dimensional Ellipse [Cascade visualizations](@ref Cascade)
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
classdef CascadeEllipse < pm.vis.Cascade

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.CascadeEllipse](@ref CascadeEllipse).<br>
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class [pm.vis.PlotEllipse](@ref PlotEllipse).<br>
        %>  \param[in]  center      :   See the corresponding input argument to the class [pm.vis.PlotEllipse](@ref PlotEllipse).<br>
        %>  \param[in]  cval        :   See the corresponding input argument to the class [pm.vis.PlotEllipse](@ref PlotEllipse).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.CascadeEllipse](@ref CascadeEllipse).<br>
        %>
        %>  \interface{CascadeEllipse}
        %>  \code{.m}
        %>
        %>      p = pm.vis.CascadeEllipse();
        %>      p = pm.vis.CascadeEllipse(gramian);
        %>      p = pm.vis.CascadeEllipse(gramian, center);
        %>      p = pm.vis.CascadeEllipse(gramian, center, cval);
        %>      p = pm.vis.CascadeEllipse(gramian, center, cval, varargin);
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
        %>  \example{CascadeEllipse}
        %>  \include{lineno} example/vis/CascadeEllipse/main.m
        %>  \vis{CascadeEllipse}
        %>  \image html example/vis/CascadeEllipse/CascadeEllipse.window.1.png width=700
        %>  \image html example/vis/CascadeEllipse/CascadeEllipse.window.2.png width=700
        %>  \image html example/vis/CascadeEllipse/CascadeEllipse.window.3.png width=700
        %>
        %>  \final{CascadeEllipse}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:18 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = CascadeEllipse(gramian, center, cval, varargin)
            %%%% Define the missing optional values as empty with the right rank.
            if  nargin < 3
                cval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end
            self = self@pm.vis.Cascade(pm.vis.PlotEllipse(gramian, center, cval), varargin{:});
        end

    end

end