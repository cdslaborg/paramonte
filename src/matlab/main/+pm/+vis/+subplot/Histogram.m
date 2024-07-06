%>  \brief
%>  This is the Histogram class for generating
%>  instances of 2-dimensional Histogram plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Histogram < pm.vis.subplot.Subplot
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
        %>  An object of [pm.vis.subplot.Histogram](@ref Histogram) class.
        %>
        %>  \interface{Histogram}
        %>  \code{.m}
        %>      s = pm.vis.subplot.Histogram(dfref);
        %>      s = pm.vis.subplot.Histogram(dfref, varargin);
        %>  \endcode
        %>
        %>
        %>  \final{Histogram}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:56 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Histogram(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Histogram", dfref, varargin{:});
        end
    end
end