%>  \brief
%>  This is the Line class for generating
%>  instances of 2-dimensional Line plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Line < pm.vis.subplot.Subplot
    methods(Access = public)
        %>
        %>  \param[in]  dfref   :   See the documentation of the corresponding input
        %>                          argument of the superclass [pm.vis.subplot.Subplot](@ref Subplot).
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.subplot.Subplot](@ref Subplot).
        %>
        %>  \return
        %>  An object of [pm.vis.subplot.Line](@ref Line) class.
        %>
        %>  \interface{Line}
        %>  \code{.m}
        %>      s = pm.vis.subplot.Line(dfref);
        %>      s = pm.vis.subplot.Line(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Line}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:59 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Line(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Line", dfref, varargin{:});
        end
    end
end