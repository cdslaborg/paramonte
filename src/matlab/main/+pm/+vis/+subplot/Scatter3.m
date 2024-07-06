%>  \brief
%>  This is the Scatter3 class for generating
%>  instances of 3-dimensional Scatter3 plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Scatter3 < pm.vis.subplot.Subplot
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
        %>  An object of [pm.vis.subplot.Scatter3](@ref Scatter3) class.
        %>
        %>  \interface{Scatter3}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Scatter3(dfref);
        %>      s = pm.vis.subplot.Scatter3(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Scatter3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Scatter3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Scatter3", dfref, varargin{:});
        end
    end
end