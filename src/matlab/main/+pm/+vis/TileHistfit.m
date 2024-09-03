%>  \brief
%>  This is the TileHistfit class for generating
%>  instances of 1-dimensional Histfit [Tile visualizations](@ref Tile)
%>  based on the relevant MATLAB
%>  intrinsic functions.
%>
%>  \note
%>  See the documentation of the attributes
%>  of the superclass [pm.vis.Tile](@ref Tile).
%>
%>  \note
%>  The input ``varargin`` can also contain the components
%>  of the ``subplot`` component of the parent object.
%>
classdef TileHistfit < pm.vis.Tile

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.TileHistfit](@ref TileHistfit).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.TileHistfit](@ref TileHistfit).<br>
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the class [pm.vis.Subplot](@ref Subplot).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \return
        %>  An object of [pm.vis.TileHistfit](@ref TileHistfit) class.
        %>
        %>  \interface{TileHistfit}
        %>  \code{.m}
        %>
        %>      t = pm.vis.TileHistfit(dfref);
        %>      t = pm.vis.TileHistfit(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{TileHistfit}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:22 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = TileHistfit(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Tile(pm.vis.SubplotHistfit(dfref), varargin{:});
        end

    end

end