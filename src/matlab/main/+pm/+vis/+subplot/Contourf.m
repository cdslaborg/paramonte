%>  \brief
%>  This is the Contourf class for generating
%>  instances of 2-dimensional Contourf plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Contourf < pm.vis.subplot.Subplot
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
        %>  An object of ``pm.vis.subplot.Contourf`` class.
        %>
        %>  \interface{Contourf}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Contourf(dfref);
        %>      s = pm.vis.subplot.Contourf(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Contourf}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:17 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Contourf(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Contourf", dfref, varargin{:});
        end
    end
end