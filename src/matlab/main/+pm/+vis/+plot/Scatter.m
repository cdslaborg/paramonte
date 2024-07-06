%>  \brief
%>  This is the Scatter class for generating
%>  instances of 2-dimensional Scatter plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Scatter < pm.vis.plot.Plot
    methods(Access = public)
        %>
        %>  \param[in]  dfref       :   See the documentation of the corresponding input
        %>                              argument of the superclass [pm.vis.plot.Plot](@ref Plot).
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \note
        %>  See the documentation of the attributes
        %>  of the superclass [pm.vis.plot.Plot](@ref Plot).
        %>
        %>  \return
        %>  An object of ``Scatter`` class.
        %>
        %>  \interface{Scatter}
        %>  \code{.m}
        %>
        %>      s = pm.vis.plot.Scatter(dfref);
        %>      s = pm.vis.plot.Scatter(dfref, varargin);
        %>
        %>  \endcode
        %>
        %>  \final{Scatter}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:41 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Scatter(dfref, varargin)
            if nargin < 1
                dfref = [];
            end
            self = self@pm.vis.plot.Plot(pm.vis.subplot.Scatter(dfref), varargin{:});
        end
    end
end