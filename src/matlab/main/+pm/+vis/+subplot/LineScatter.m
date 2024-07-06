%>  \brief
%>  This is the LineScatter class for generating
%>  instances of 2-dimensional LineScatter plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef LineScatter < pm.vis.subplot.Subplot
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
        %>  An object of [pm.vis.subplot.LineScatter](@ref LineScatter) class.
        %>
        %>  \interface{LineScatter}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.LineScatter(dfref);
        %>      s = pm.vis.subplot.LineScatter(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{LineScatter}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 6:00 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = LineScatter(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("LineScatter", dfref, varargin{:});
        end
    end
end