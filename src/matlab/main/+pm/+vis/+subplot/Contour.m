%>  \brief
%>  This is the Contour class for generating
%>  instances of 2-dimensional Contour plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Contour < pm.vis.subplot.Subplot
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
        %>  An object of [pm.vis.subplot.Contour](@ref Contour) class.
        %>
        %>  \interface{Contour}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Contour(dfref);
        %>      s = pm.vis.subplot.Contour(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Contour}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:16 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Contour(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Contour", dfref, varargin{:});
        end
    end
end