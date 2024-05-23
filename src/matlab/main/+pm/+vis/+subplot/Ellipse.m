%>  \brief
%>  This is the Ellipse class for generating
%>  instances of 2-dimensional Ellipse plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Ellipse < pm.vis.subplot.Ellipse3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>
        %>  \details
        %>  In the following documentation,
        %>  <ol>
        %>      <li>  the variable ``ndim`` represents the number of ellipsoids in the input data.
        %>          The value of ``ndim`` is inferred from the shapes
        %>          of the input ``gramian`` and ``center`` arguments
        %>          as ``max([2, size(gramian, 1), size(center, 1)])``.
        %>
        %>      <li>  the variable ``nell`` represents the number of ellipsoids in the input data.
        %>          The value of ``nell`` is inferred from the shapes
        %>          of the input ``gramian``, ``center``, and ``cval`` arguments
        %>          as ``max([size(gramian, 3), size(center, 2), size(cval, 2)])``.
        %>          If the above expression yields zero, then ``nell`` is set to ``75``.
        %>
        %>      <li>  the variable ``npnt`` represents the number of points used in visualizing each ellipsoid.
        %>          The value of ``npnt`` is inferred from the shapes of the input argument ``cval`` as ``size(cval, 1)``.
        %>          If the above expression yields zero or one, then ``npnt`` is set to ``100``.
        %>
        %>  </ol>
        %>
        %>  \note
        %>  This class is merely a simple extension of
        %>  the superclass ``pm.vis.subplot.Ellipse3``.
        %>  The only difference is that all ``zval`` values for
        %>  all points of all coordinates are set to default values.
        %>  The only tangible action this class performs is
        %>  to adjust the camera view to a 2D axis.
        %>
        %>  \param[in]  gramian     :   See the corresponding input argument to the constructor of the parent object.
        %>  
        %>  \param[in]  center      :   See the corresponding input argument to the constructor of the parent object.
        %>  
        %>  \param[in]  cval        :   See the corresponding input argument to the constructor of the parent object.
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  See below and also the documentation of the attributes
        %>  of the superclass ``pm.vis.subplot.Ellipse3``.
        %>
        %>  \return
        %>  An object of class ``pm.vis.subplot.Ellipse``.
        %>
        %>  \interface{Ellipse}
        %>  \code{.m}
        %>
        %>      s = pm.vis.subplot.Ellipse();
        %>      s = pm.vis.subplot.Ellipse(gramian);
        %>      s = pm.vis.subplot.Ellipse(gramian, center);
        %>      s = pm.vis.subplot.Ellipse(gramian, center, cval);
        %>      s = pm.vis.subplot.Ellipse(gramian, center, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Ellipse}
        %>
        %>      s = pm.vis.subplot.Ellipse();
        %>      s.make("dims", [1, 2]);
        %>
        %>  \final{Ellipse}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:21 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Ellipse(gramian, center, cval, varargin)

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
            self = self@pm.vis.subplot.Ellipse3(gramian, center, [], cval, varargin{:})

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Generate a plot from the selected dimensions
        %>  of the input ellipsoid data to the object constructor.
        %>
        %>  \note
        %>  This method is pure, except for the changes in ``fout`` component.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>  \endcode
        %>
        %>      pm.vis.subplot.Ellipse.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>
        %>      s = pm.vis.subplot.Ellipse();
        %>      s.make("dims", [1, 2]);
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 10:23 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            make@pm.vis.subplot.Ellipse3(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            view(2);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end