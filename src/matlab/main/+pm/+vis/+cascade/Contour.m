%>  \brief
%>  This is the Contour class for generating
%>  instances of 2-dimensional Contour plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Contour < pm.vis.cascade.Cascade
    methods(Access = public)
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the class ``pm.vis.plot.Plot``.
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``template`` component of the parent object.
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass ``pm.vis.cascade.Cascade``.
        %>
        %>  \return
        %>  An object of ``pm.vis.cascade.Contour`` class.
        %>
        %>  \interface{Contour}
        %>  \code{.m}
        %>
        %>      p = pm.vis.cascade.Contour(dfref);
        %>      p = pm.vis.cascade.Contour(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Contour}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:12 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Contour(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.cascade.Cascade(pm.vis.plot.Contour(dfref), varargin{:});
        end
    end
end