%>  \brief
%>  This is the LineScatter class for generating
%>  instances of 2-dimensional LineScatter tiles
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef LineScatter < pm.vis.tile.Tile
    methods(Access = public)
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the superclass [pm.vis.tile.Tile](@ref Tile).
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
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.tile.Tile](@ref Tile).
        %>
        %>  \return
        %>  An object of [pm.vis.tile.LineScatter](@ref LineScatter) class.
        %>
        %>  \interface{LineScatter}
        %>  \code{.m}
        %>
        %>      t = pm.vis.tile.LineScatter(dfref);
        %>      t = pm.vis.tile.LineScatter(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{LineScatter}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:39 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = LineScatter(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.tile.Tile(pm.vis.subplot.LineScatter(dfref), varargin{:});
        end
    end
end