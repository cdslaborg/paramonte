%>  \brief
%>  This is the TileLine3 class for generating
%>  instances of 3-dimensional Line3 [Tile visualizations](@ref Tile)
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef TileLine3 < pm.vis.Tile

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.TileLine3](@ref TileLine3).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.TileLine3](@ref TileLine3).<br>
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the class [pm.vis.Subplot](@ref Subplot).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
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
        %>  of the superclass [pm.vis.Tile](@ref Tile).
        %>
        %>  \return
        %>  An object of [pm.vis.TileLine3](@ref TileLine3) class.
        %>
        %>  \interface{TileLine3}
        %>  \code{.m}
        %>      t = pm.vis.TileLine3(dfref);
        %>      t = pm.vis.TileLine3(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{TileLine3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:29 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = TileLine3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.Tile(pm.vis.SubplotLine3(dfref), varargin{:});
        end

    end

end