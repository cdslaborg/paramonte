%>  \brief
%>  This is the Ellipse3 class for generating
%>  instances of 3-dimensional Ellipse3 plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Ellipse3 < pm.vis.cascade.Cascade
    methods(Access = public)
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the class ``pm.vis.plot.Ellipse3``.
        %>  
        %>  \param[in]  center      :   See the corresponding input argument to the class ``pm.vis.plot.Ellipse3``.
        %>  
        %>  \param[in]  zval        :   See the corresponding input argument to the class ``pm.vis.plot.Ellipse3``.
        %>  
        %>  \param[in]  cval        :   See the corresponding input argument to the class ``pm.vis.plot.Ellipse3``.
        %>  
        %>  \param[in]  varargin    :    Any ``property, value`` pair of the parent object.
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
        %>  See below and also the documentation of the
        %>  attributes of the superclass ``pm.vis.figure.Figure``.
        %>
        %>  \return
        %>  An object of class ``pm.vis.cascade.Ellipse3``.
        %>
        %>  \interface{Ellipse3}
        %>  \code{.m}
        %>
        %>      p = pm.vis.cascade.Ellipse3();
        %>      p = pm.vis.cascade.Ellipse3(gramian);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center, zval);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center, zval, cval);
        %>      p = pm.vis.cascade.Ellipse3(gramian, center, zval, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Ellipse3}
        %>
        %>      p = pm.vis.cascade.Ellipse3();
        %>      p.make("dims", [1, 2]);
        %>
        %>  \final{Ellipse3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 7:21 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Ellipse3(gramian, center, zval, cval, varargin)
            %%%% Define the missing optional values as empty with the right rank.
            if  nargin < 4
                cval = zeros(0, 0);
            end
            if  nargin < 3
                zval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end
            self = self@pm.vis.cascade.Cascade(pm.vis.plot.Ellipse3(gramian, center, zval, cval), varargin{:});
        end
    end
end