%>  \brief
%>  This is the Histfit class for generating
%>  instances of 2-dimensional Histfit plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Histfit < pm.vis.subplot.Subplot
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
        %>  An object of [pm.vis.subplot.Histfit](@ref Histfit) class.
        %>
        %>  \interface{Histfit}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Histfit(dfref);
        %>      s = pm.vis.subplot.Histfit(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Histfit}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:55 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Histfit(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Histfit", dfref, varargin{:});
        end
    end
end