%>  \brief
%>  This is the Ellipse class for generating
%>  instances of 2-dimensional Ellipse tiles
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Ellipse < pm.vis.tile.Tile
    methods(Access = public)
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class ``pm.vis.subplot.Ellipse``.
        %>  
        %>  \param[in]  center      :   See the corresponding input argument to the class ``pm.vis.subplot.Ellipse``.
        %>  
        %>  \param[in]  cval        :   See the corresponding input argument to the class ``pm.vis.subplot.Ellipse``.
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
        %>  See below and also the documentation of the
        %>  attributes of the superclass ``pm.vis.tile.Tile``.
        %>
        %>  \return
        %>  An object of class ``pm.vis.tile.Ellipse``.
        %>
        %>  \interface{Ellipse}
        %>  \code{.m}
        %>  \endcode
        %>
        %>      t = pm.vis.tile.Ellipse();
        %>      t = pm.vis.tile.Ellipse(gramian);
        %>      t = pm.vis.tile.Ellipse(gramian, center);
        %>      t = pm.vis.tile.Ellipse(gramian, center, cval);
        %>      t = pm.vis.tile.Ellipse(gramian, center, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Ellipse}
        %>
        %>      t = pm.vis.tile.Ellipse();
        %>      t.make("dims", [1, 2]);
        %>
        %>  \final{Ellipse}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:16 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Ellipse(gramian, center, cval, varargin)
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
            self = self@pm.vis.tile.Tile(pm.vis.subplot.Ellipse(gramian, center, cval), varargin{:});
        end
    end
end