%>  \brief
%>  This is the SubplotEllipse class for generating
%>  instances of 2-dimensional Ellipse [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB
%>  intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotEllipse](@ref SubplotEllipse::SubplotEllipse) for example usage.<br>
%>
%>  \note
%>  See the documentation of the attributes
%>  of the superclass [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
%>
%>  \see
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Corner](@ref Corner)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef SubplotEllipse < pm.vis.SubplotEllipse3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotEllipse](@ref SubplotEllipse).<br>
        %>
        %>  \details
        %>  In the following documentation,
        %>  <ol>
        %>      <li>    The variable ``ndim`` represents the number of ellipsoids in the input data.<br>
        %>              The value of ``ndim`` is inferred from the shapes of the input ``gramian`` and
        %>              ``center`` arguments as ``max([2, size(gramian, 1), size(center, 1)])``.<br>
        %>
        %>      <li>    The variable ``nell`` represents the number of ellipsoids in the input data.<br>
        %>              The value of ``nell`` is inferred from the shapes of the input ``gramian``, ``center``,
        %>              and ``cval`` arguments as ``max([size(gramian, 3), size(center, 2), size(cval, 2)])``.<br>
        %>              If the above expression yields zero, then ``nell`` is set to ``75``.<br>
        %>
        %>      <li>    The variable ``npnt`` represents the number of points used in visualizing each ellipsoid.<br>
        %>              The value of ``npnt`` is inferred from the shapes of the input argument ``cval`` as ``size(cval, 1)``.<br>
        %>              If the above expression yields zero or one, then ``npnt`` is set to ``100``.<br>
        %>
        %>  </ol>
        %>
        %>  \note
        %>  This class is merely a simple extension of the superclass [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>  The only difference is that all ``zval`` values for all points of all coordinates are set to default values.<br>
        %>  The only tangible action this class performs is to adjust the camera view to a 2D axis.<br>
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the constructor of the parent object.<br>
        %>  \param[in]  center      :   See the corresponding input argument to the constructor of the parent object.<br>
        %>  \param[in]  cval        :   See the corresponding input argument to the constructor of the parent object.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output object of class [pm.vis.SubplotEllipse](@ref SubplotEllipse).<br>
        %>
        %>  \interface{SubplotEllipse}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotEllipse();
        %>      s = pm.vis.SubplotEllipse(gramian);
        %>      s = pm.vis.SubplotEllipse(gramian, center);
        %>      s = pm.vis.SubplotEllipse(gramian, center, cval);
        %>      s = pm.vis.SubplotEllipse(gramian, center, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>
        %>  \example{SubplotEllipse}
        %>  \include{lineno} example/vis/SubplotEllipse/main.m
        %>  \vis{SubplotEllipse}
        %>  <br>\image html example/vis/SubplotEllipse/SubplotEllipse.1.png width=700
        %>
        %>  \final{SubplotEllipse}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:21 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotEllipse(gramian, center, cval, varargin)

            %%%% Define the missing optional values as empty with the right rank.

            if  nargin < 3
                cval = zeros(0, 0);
            end
            if  nargin < 2
                center = zeros(0, 0);
            end
            if  nargin < 1
                gramian = zeros(0, 0, 0);
            end
            self = self@pm.vis.SubplotEllipse3(gramian, center, [], cval, varargin{:})

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Generate a plot from the selected dimensions
        %>  of the input ellipsoid data to the object constructor.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.SubplotEllipse](@ref SubplotEllipse)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>  \endcode
        %>
        %>      pm.vis.SubplotEllipse.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>  \include{lineno} example/vis/SubplotEllipse/main.m
        %>  \vis{make}
        %>  <br>\image html example/vis/SubplotEllipse/SubplotEllipse.1.png width=700
        %>
        %>  \final{make}
        %>
        %>  \note
        %>  This method is pure, except for the changes in ``fout`` component.<br>
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:23 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            make@pm.vis.SubplotEllipse3(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            view(2);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end