%>  \brief
%>  This is the [pm.vis.TileEllipse](@ref TileEllipse) class for generating
%>  instances of 2-dimensional Ellipse [Tile visualizations](@ref Tile)
%>  based on the relevant MATLAB intrinsic functions.<br>
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
classdef TileEllipse < pm.vis.Tile

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.TileEllipse](@ref TileEllipse).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.TileEllipse](@ref TileEllipse).<br>
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class [pm.vis.SubplotEllipse](@ref SubplotEllipse).
        %>  \param[in]  center      :   See the corresponding input argument to the class [pm.vis.SubplotEllipse](@ref SubplotEllipse).
        %>  \param[in]  cval        :   See the corresponding input argument to the class [pm.vis.SubplotEllipse](@ref SubplotEllipse).
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   An object of [pm.vis.TileEllipse](@ref TileEllipse) class.<br>
        %>
        %>  \interface{TileEllipse}
        %>  \code{.m}
        %>
        %>      t = pm.vis.TileEllipse();
        %>      t = pm.vis.TileEllipse(gramian);
        %>      t = pm.vis.TileEllipse(gramian, center);
        %>      t = pm.vis.TileEllipse(gramian, center, cval);
        %>      t = pm.vis.TileEllipse(gramian, center, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.<br>
        %>
        %>  \note
        %>  See below and also the documentation of the
        %>  attributes of the superclass [pm.vis.Tile](@ref Tile).<br>
        %>
        %>  \example{TileEllipse}
        %>
        %>      t = pm.vis.TileEllipse();
        %>      t.make("dimx", 1, "dimy", 2);
        %>
        %>  \final{TileEllipse}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:16 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = TileEllipse(gramian, center, cval, varargin)
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
            self = self@pm.vis.Tile(pm.vis.SubplotEllipse(gramian, center, cval), varargin{:});
        end
    end
end