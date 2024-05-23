%>  \brief
%>  This is the Line3 class for generating
%>  instances of 3-dimensional Line3 plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Line3 < pm.vis.subplot.Subplot
    methods(Access = public)
        %>
        %>  \param[in]  dfref   :   See the documentation of the corresponding input
        %>                          argument of the superclass ``pm.vis.subplot.Subplot``.
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass ``pm.vis.subplot.Subplot``.
        %>
        %>  \return
        %>  An object of ``pm.vis.subplot.Line3`` class.
        %>
        %>  \interface{Line3}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Line3(dfref);
        %>      s = pm.vis.subplot.Line3(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Line3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:01 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Line3(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.subplot.Subplot("Line3", dfref, varargin{:});
        end
    end
end