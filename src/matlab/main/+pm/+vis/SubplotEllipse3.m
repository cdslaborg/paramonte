%>  \brief
%>  This is the SubplotEllipse3 class for generating
%>  instances of 3-dimensional Ellipse [Subplot visualizations](@ref Subplot)
%>  based on the relevant MATLAB intrinsic functions.<br>
%>
%>  \note
%>  See the documentation of the constructor of the class
%>  [pm.vis.SubplotEllipse3](@ref SubplotEllipse3::SubplotEllipse3) for example usage.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3).<br>
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
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef SubplotEllipse3 < pm.vis.SubplotLineScatter3
    properties(Access = public)
        %>
        %>  ``names``
        %>
        %>  The vector of MATLAB strings containing the names of dimensions
        %>  of the space within which the Gramian matrices are defined.<br>
        %>  The default is ``"Dimension i"`` where ``i`` is replaced
        %>  by the ID of the corresponding axis.<br>
        %>
        names = [];
        %>
        %>  ``dimx``
        %>
        %>  The MATLAB scalar (or vector) of whole-number(s) representing the dimension(s) of the
        %>  input ``gramian`` and ``center`` data to be used on the x-axis of 2D ellipsoid visualization.<br>
        %>  The default is ``1`` if ``dimy`` component is also unspecified, or the appropriate data column adjacent to ``dimmy``.<br>
        %>  If the condition ``dimx == dimy`` holds, then the resulting ellipses will be circles.<br>
        %>
        dimx = [];
        %>
        %>  ``dimy``
        %>
        %>  The MATLAB scalar (or vector) of whole-number(s) representing the dimension(s) of the
        %>  input ``gramian`` and ``center`` data to be used on the y-axis of 2D ellipsoid visualization.<br>
        %>  The default is ``1`` if ``dimx`` component is also unspecified, or the appropriate data column adjacent to ``dimmx``.<br>
        %>  If the condition ``dimx == dimy`` holds, then the resulting ellipses will be circles.<br>
        %>
        dimy = [];
        %>
        %>  ``ellindex``
        %>
        %>  The vector of MATLAB whole numbers representing the ellindex of the 2D
        %>  ellipsoids to display, represented by the input ``gramian`` and ``center``.<br>
        %>  The specified ellindex will serve as the visualization values on the z-axis.<br>
        %>  The default is ``ellindex = pm.array.logrange(start = 1, stop = nell, count = 75);``,
        %>  where ``nell`` represents the number of ellipsoids to visualize.<br>
        %>
        ellindex = [];
        %>
        %>  ``gramian``
        %>
        %>  A scalar object of class [pm.data.DataRef](@ref DataRef)
        %>  containing the user-specified Gramian data to visualize.<br>
        %>
        gramian = [];
        %>
        %>  ``center``
        %>
        %>  A scalar object of class [pm.data.DataRef](@ref DataRef)
        %>  containing the user-specified center data to visualize.<br>
        %>
        center = [];
        %>
        %>  ``zval``
        %>
        %>  A scalar object of class [pm.data.DataRef](@ref DataRef)
        %>  containing the user-specified z-axis data to visualize.<br>
        %>
        zval = [];
        %>
        %>  ``cval``
        %>
        %>  A scalar object of class [pm.data.DataRef](@ref DataRef)
        %>  containing the user-specified color data in colormap.<br>
        %>  If empty or unspecified, it will be set to ``zval`` component.<br>
        %>
        cval = [];
    end

    properties(Access = public, Hidden)
        %>
        %>  ``ndim``
        %>
        %>  A scalar MATLAB whole number representing the total number of the
        %>  dimensions of ellipsoids identified based on the input that can be visualized.<br>
        %>  It is set to ``max(2, size(gramian, 2), size(center, 1))``.<br>
        %>
        ndim = [];
        %>
        %>  ``nell``
        %>
        %>  A scalar MATLAB whole number representing the total number of
        %>  ellipsoids identified based on the input that can be visualized.<br>
        %>  It is set to ``max(20, size(gramian, 3), size(center, 2), size(zval, 2))``.<br>
        %>
        nell = [];
        %>
        %>  ``npnt``
        %>
        %>  The scalar MATLAB whole number representing the number of
        %>  points with which the 2D ellipsoid boundaries are delineated.<br>
        %>  It is set to ``max(size(zval, 1), size(cval, 1))`` or otherwise, ``100``.<br>
        %>
        npnt = [];
        %>
        %>  ``workspace``
        %>
        %>  A scalar MATLAB ``struct()`` containing the current visualization information.<br>
        %>  This information is updated with every call to the ``make()`` method.<br>
        %>  The contents of this component are temporary and must remain read-only.<br>
        %>  This component is solely provided for better insight into the internal
        %>  workings of the ``make()`` method of the parent object and
        %>  the resulting visualization.<br>
        %>
        workspace = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.vis.SubplotEllipse](@ref SubplotEllipse).<br>
        %>
        %>  \details
        %>  In the following documentation,
        %>  <ol>
        %>      <li>    The variable ``ndim`` represents the number of ellipsoids in the input data.<br>
        %>              The value of ``ndim`` is inferred from the shapes of the input ``gramian``
        %>              and ``center`` arguments as ``max([2, size(gramian, 1), size(center, 1)])``.<br>
        %>
        %>      <li>    The variable ``nell`` represents the number of ellipsoids in the input data.<br>
        %>              The value of ``nell`` is inferred from the shapes
        %>              of the input ``gramian``, ``center``, ``zval``, and ``cval`` arguments
        %>              as ``max([size(gramian, 3), size(center, 2), size(zval, 2), size(cval, 2)])``.<br>
        %>              If the above expression yields zero, then ``nell`` is set to ``75``.<br>
        %>
        %>      <li>    The variable ``npnt`` represents the number of points used in visualizing each ellipsoid.<br>
        %>              The value of ``npnt`` is inferred from the shapes of the input arguments
        %>              ``zval`` and ``czval`` as ``max(size(zval, 1), size(cval, 1))``.<br>
        %>              If the above expression yields zero or one, then ``npnt`` is set to ``100``.<br>
        %>  </ol>
        %>
        %>  \param[in]  gramian     :   The MATLAB (function handle returning an) object that can be either,<br>
        %>                              <ol>
        %>                                  <li>    an empty object (e.g., ``[]``), in which case,
        %>                                          the value ``repmat(eye(ndim, ndim), 1, 1, nell)`` will be used.<br>
        %>
        %>                                  <li>    a scalar, containing the value of the diagonal elements
        %>                                          of the diagonal matrices representing the Gramian matrices
        %>                                          of all ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a vector of shape ``[1, 1, nell]`` such that ``gramian(1, 1, iell)``
        %>                                          represents the diagonal value to be used in constructing the Gramian
        %>                                          matrix of the ellipsoid ``iell` in the input ``nell`` ellipsoid data.<br>
        %>
        %>                                  <li>    a matrix of shape ``[ndim, ndim, 1]`` used in constructing the Gramian
        %>                                          matrices of all ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a tensor of shape ``[ndim, ndim, nell]`` such that ``gramian(:, :, iell)``
        %>                                          represents the Gramian matrix of the ellipsoid ``iell` in the input data.<br>
        %>                              </ol>
        %>                              While it is possible to pass a data directly to this class,
        %>                              it is highly recommended to pass a function handle that returns
        %>                              such data when called, allowing the visualization data to be
        %>                              dynamically updated when needed.<br>
        %>  \param[in]  center      :   The MATLAB (function handle returning an) object that can be either,<br>
        %>                              <ol>
        %>                                  <li>    an empty object (e.g., ``[]``), in which case,
        %>                                          the value ``zeros(ndim, nell)`` will be used.<br>
        %>
        %>                                  <li>    a scalar, containing the value to assign to all coordinates
        %>                                          of the centers of all ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a vector of shape ``[1, nell]`` such that ``center(1, iell)``
        %>                                          represents the value of all coordinates of the center of the
        %>                                          ellipsoid ``iell` in the input ``nell`` ellipsoid data.<br>
        %>
        %>                                  <li>    a vector of shape ``[ndim, 1]`` such that ``center(idim, 1)``
        %>                                          represents the value of ``idim`` coordinate of the centers of
        %>                                          all ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a matrix of shape ``[ndim, nell]`` such that ``center(:, iell)`` represents
        %>                                          the coordinates of the center of the ellipsoid ``iell` in the input data.<br>
        %>                              </ol>
        %>                              While it is possible to pass data directly to this class,
        %>                              it is highly recommended to pass a function handle that returns
        %>                              such data when called, allowing the visualization data to be
        %>                              dynamically updated when needed.<br>
        %>  \param[in]  zval        :   The MATLAB (function handle returning an) object that can be either,<br>
        %>                              <ol>
        %>                                  <li>    an empty object (e.g., ``[]``), in which case,
        %>                                          the value ``repmat(1 : nell, npnt, 1)`` will be used.<br>
        %>
        %>                                  <li>    a scalar, containing the value to assign to the z-axis
        %>                                          coordinates of all ``npnt`` points used in the visualization
        %>                                          of ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a vector of shape ``[1, nell]`` such that ``zval(1, iell)``
        %>                                          represents the value of z-axis coordinates of all ``npnt``
        %>                                          points used in the visualization of the ellipsoid
        %>                                          ``iell` in the input ``nell`` ellipsoid data.<br>
        %>
        %>                                  <li>    a vector of shape ``[npnt, 1]`` such that ``zval(ipnt, 1)``
        %>                                          contains the z-axis value of ``ipnt`` points in the representations
        %>                                          of all ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a matrix of shape ``[npnt, nell]`` such that ``zval(:, iell)`` represents
        %>                                          the z-axis coordinates of all ``npnt`` points used in the visualization of
        %>                                          the ellipsoid ``iell` in the input data.<br>
        %>                              </ol>
        %>                              While it is possible to pass data directly to this class,
        %>                              it is highly recommended to pass a function handle that returns
        %>                              such data when called, allowing the visualization data to be
        %>                              dynamically updated when needed.<br>
        %>  \param[in]  cval        :   The MATLAB (function handle returning an) object that can be either,<br>
        %>                              <ol>
        %>                                  <li>    an empty object (e.g., ``[]``), in which case,
        %>                                          the input value ``zval`` will be used.<br>
        %>
        %>                                  <li>    a scalar, containing the value to assign to color values
        %>                                          of all ``npnt`` points used in the visualization
        %>                                          of  ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a vector of shape ``[1, nell]`` such that ``cval(1, iell)``
        %>                                          represents the color values of all ``npnt`` points used in the
        %>                                          visualization of the ellipsoid ``iell` in the input ``nell`` ellipsoid data.<br>
        %>
        %>                                  <li>    a vector of shape ``[npnt, 1]`` such that ``cval(ipnt, 1)``
        %>                                          contains the color values of ``ipnt`` points in the representations
        %>                                          of all ``nell`` ellipsoids in the input data.<br>
        %>
        %>                                  <li>    a matrix of shape ``[npnt, nell]`` such that ``center(:, iell)`` represents
        %>                                          the color values of all ``npnt`` points used in the visualization of
        %>                                          the ellipsoid ``iell` in the input data.<br>
        %>                              </ol>
        %>                              While it is possible to pass data directly to this class,
        %>                              it is highly recommended to pass a function handle that returns
        %>                              such data when called, allowing the visualization data to be
        %>                              dynamically updated when needed.<br>
        %>                              The input value for ``cval`` is used only if ``colormap.enabled`` component of the
        %>                              output object of class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3) is set to ``true``.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output object of class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3).<br>
        %>
        %>  \interface{SubplotEllipse3}
        %>  \code{.m}
        %>      s = pm.vis.SubplotEllipse3();
        %>      s = pm.vis.SubplotEllipse3(gramian);
        %>      s = pm.vis.SubplotEllipse3(gramian, center);
        %>      s = pm.vis.SubplotEllipse3(gramian, center, zval);
        %>      s = pm.vis.SubplotEllipse3(gramian, center, zval, cval);
        %>      s = pm.vis.SubplotEllipse3(gramian, center, zval, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  To generate a 2D plot of ellipsoids, simply execute the
        %>  MATLAB command ``view(2)`` to change the camera view to 2D.<br>
        %>
        %>  \note
        %>  See also the documentation of the attributes
        %>  of the superclass [pm.vis.SubplotLineScatter3](@ref SubplotLineScatter3).<br>
        %>
        %>  \example{SubplotEllipse3}
        %>  \include{lineno} example/vis/SubplotEllipse3/main.m
        %>  \vis{SubplotEllipse3}
        %>  <br>\image html example/vis/SubplotEllipse3/SubplotEllipse3.1.png width=700
        %>
        %>  \final{SubplotEllipse3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:35 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SubplotEllipse3(gramian, center, zval, cval, varargin)

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

            varargin =  { "gramian", pm.data.DataRef(gramian) ...
                        , "center", pm.data.DataRef(center) ...
                        , "zval", pm.data.DataRef(zval) ...
                        , "cval", pm.data.DataRef(cval) ...
                        , varargin{:} ...
                        };
            self = self@pm.vis.SubplotLineScatter3([], varargin{:})

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Reset the properties of the plot to the original default settings.<br>
        %>
        %>  \details
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.SubplotEllipse3.reset() % reset the plot to the default settings.
        %>      pm.vis.SubplotEllipse3.reset(varargin)
        %>
        %>  \endcode
        %>
        %>  \final{len}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 5:42 PM, University of Texas at Arlington<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)
            self.dimx = [];
            self.dimy = [];
            self.npnt = [];
            self.names = [];
            self.ellindex = [];
            reset@pm.vis.SubplotLineScatter3(self, varargin{:});
            self.surface.lineWidth = [];
            self.scatter3.enabled = [];
            self.surface.enabled = [];
            self.plot3.enabled = [];
            self.colormap.map = []; %"winter"
            self.colormap.enabled = [];
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %>  \brief
        %>  Prepare the subplot for visualization.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``premake()`` method.<br>
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      pm.vis.SubplotEllipse3.premake();
        %>      pm.vis.SubplotEllipse3.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{premake}
        %>  \code{.m}
        %>
        %>      s = pm.vis.SubplotEllipse3();
        %>      s.premake("dimx", 1, "dimy", 2], "colormap", {"map", "autumn"})
        %>
        %>  \endcode
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 5:45 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function premake(self, varargin)

            if ~isempty(varargin)
                self.hash2comp(varargin); % parse arguments
            end

            %%%%
            %%%% Set the visualization types.
            %%%%

            if  isempty(self.colormap.map)
                self.colormap.map = "winter";
            end

            if  isempty(self.colormap.enabled)
                self.colormap.enabled = true;
            end

            if  isempty(self.colorbar.enabled)
                self.colorbar.enabled = self.colormap.enabled;
            end

            if  isempty(self.plot3.enabled)
                self.plot3.enabled = ~self.colormap.enabled;
            end

            if  isempty(self.surface.enabled)
                self.surface.enabled = self.colormap.enabled;
            end

            if  isempty(self.scatter3.enabled)
                self.scatter3.enabled = ~(self.plot3.enabled || self.surface.enabled);
            end

            %%%%
            %%%% Set the visualization specs.
            %%%%

            if  isempty(self.surface.lineWidth)
                self.surface.lineWidth = 1.5;
            end

            %%%%
            %%%% Set the visualization specs of the parent (and set the data components).
            %%%%

            premake@pm.vis.SubplotLineScatter3(self);

            %%%%
            %%%% Set the line color to default MATLAB blue only if scatter is disabled and colormap is disabled and color is unset by the user.
            %%%%

            if ~(self.colormap.enabled || self.scatter3.enabled)
                [val, failed] = pm.matlab.hashmap.getKeyVal("plot3", varargin);
                if ~failed
                    [~, failed] = pm.matlab.hashmap.getKeyVal("color", val);
                end
                if  failed
                    self.plot3.color = [0, 0.4470, 0.7410];
                end
            end

            %%%%
            %%%% Set the data.
            %%%%

            self.workspace = struct();
            self.workspace.cval = self.cval.copy();
            self.workspace.zval = self.zval.copy();
            self.workspace.center = self.center.copy();
            self.workspace.gramian = self.gramian.copy();

            %%%%
            %%%% Get the shape of input data.
            %%%%

            %%%% The following ``size()`` syntax is essential for compatibility with MATLAB R2019a.

            self.workspace.shape = struct();

            [s1, s2] = size(self.workspace.cval);
            self.workspace.shape.cval = [s1, s2];

            [s1, s2] = size(self.workspace.zval);
            self.workspace.shape.zval = [s1, s2];

            [s1, s2] = size(self.workspace.center);
            self.workspace.shape.center = [s1, s2];

            [s1, s2, s3] = size(self.workspace.gramian);
            self.workspace.shape.gramian = [s1, s2, s3];

            %%%%
            %%%% Infer ``ndim``.
            %%%%

            self.ndim = max([2, self.workspace.shape.gramian(1), self.workspace.shape.center(1)]);
            if  self.ndim < 2
                self.ndim = 2;
            end

            %%%%
            %%%% Infer ``nell``.
            %%%%

            self.nell = max([self.workspace.shape.gramian(3), self.workspace.shape.center(2), self.workspace.shape.zval(2), self.workspace.shape.cval(2)]);
            if  self.nell == 0
                self.nell = 75;
            end

            %%%%
            %%%% Infer ``npnt``.
            %%%%

            if  isempty(self.npnt)
                self.npnt = max(self.workspace.shape.zval(1), self.workspace.shape.cval(1));
                if  self.npnt < 2
                    self.npnt = 100;
                end
            end

            %%%%
            %%%% Assert data shape consistencies.
            %%%%

            asserted = true;
            asserted = asserted && self.ndim == self.workspace.shape.center(1) || self.workspace.shape.center(1) < 2;
            asserted = asserted && self.ndim == self.workspace.shape.gramian(1) || self.workspace.shape.gramian(1) < 2;
            asserted = asserted && self.ndim == self.workspace.shape.gramian(2) || self.workspace.shape.gramian(2) < 2;
            asserted = asserted && self.nell == self.workspace.shape.gramian(3) || self.workspace.shape.gramian(3) < 2;
            asserted = asserted && self.nell == self.workspace.shape.center(2) || self.workspace.shape.center(2) < 2;
            asserted = asserted && self.nell == self.workspace.shape.zval(2) || self.workspace.shape.zval(2) < 2;
            asserted = asserted && self.nell == self.workspace.shape.cval(2) || self.workspace.shape.cval(2) < 2;

            if ~asserted
                help("pm.vis.SubplotEllipse3")
                disp("size(gramian)")
                disp( self.workspace.shape.gramian )
                disp("size(center)")
                disp( self.workspace.shape.center )
                disp("size(zval)")
                disp( self.workspace.shape.zval )
                disp("size(cval)")
                disp( self.workspace.shape.cval )
                disp("[ndim, nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``gramian``, ``center``, ``zval``, and ``cval`` are incompatible." + newline ...
                        + "For more information, see the documentation of the input arguments to the object constructor displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Set the default dimx and dimy or verify the custom values to visualize.
            %%%%

            if  isempty(self.dimx) && isempty(self.dimy)
                self.workspace.dimx = 1;
                self.workspace.dimy = 2;
            elseif isempty(self.dimx) && ~isempty(self.dimy)
                self.workspace.dimy = self.dimy;
                self.workspace.dimx = zeros(numel(self.dimy), 1);
                for i = 1 : numel(self.dimy)
                    if  self.dimy(i) == 1
                        self.workspace.dimx(i) = self.dimy(i) + 1;
                    else
                        self.workspace.dimx(i) = self.dimy(i) - 1;
                    end
                end
            elseif isempty(self.dimy) && ~isempty(self.dimx)
                self.workspace.dimx = self.dimx;
                self.workspace.dimy = zeros(numel(self.dimx), 1);
                for i = 1 : numel(self.dimx)
                    if  self.dimx(i) == 1
                        self.workspace.dimy(i) = self.dimx(i) + 1;
                    else
                        self.workspace.dimy(i) = self.dimx(i) - 1;
                    end
                end
            else
                self.workspace.dimx = self.dimx;
                self.workspace.dimy = self.dimy;
            end

            if  numel(self.workspace.dimx) ~= numel(self.workspace.dimy)
                if  numel(self.workspace.dimx) == 1
                    self.workspace.dimx = self.workspace.dimx * ones(numel(self.workspace.dimy), 1);
                elseif numel(self.workspace.dimy) == 1
                    self.workspace.dimy = self.workspace.dimy * ones(numel(self.workspace.dimx), 1);
                else
                    help("pm.vis.SubplotEllipse3")
                    disp("self.dimx")
                    disp( self.dimx )
                    disp("size(self.dimx)")
                    disp( size(self.dimx) )
                    disp("self.dimy")
                    disp( self.dimy )
                    disp("size(self.dimy)")
                    disp( size(self.dimy) )
                    error   ( newline ...
                            + "The sizes of the specified ``dimx`` and ``dimy`` must either be equal or one must have size of ``1``." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end
            end

           %%%%
           %%%% Check the consistency of shape of ``self.names``.
           %%%%

            if ~isempty(self.names) && pm.array.len(self.names) ~= self.ndim
                help("pm.vis.SubplotEllipse3")
                disp("string(self.names)")
                disp( string(self.names) )
                disp("length(string(self.names))")
                disp( length(string(self.names)) )
                error   ( newline ...
                        + "The component ``names`` must be a vector of size ``ndim`` of strings." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Set the visualization axes names.
            %%%%

            if ~isempty(self.names)
                self.workspace.names = string(self.names);
                self.workspace.names = self.workspace.names(:);
            else
                self.workspace.names = strings(self.ndim, 1);
                for idim = 1 : self.ndim
                    self.workspace.names(idim) = "Dimension " + string(idim);
                end
            end

            if ~pm.array.len(self.xlabel.txt)
                self.xlabel.txt = string(self.workspace.names(self.workspace.dimx(:)));
            end
            if ~pm.array.len(self.ylabel.txt)
                self.ylabel.txt = string(self.workspace.names(self.workspace.dimy(:)));
            end
            if ~pm.array.len(self.zlabel.txt)
                if ~isempty(self.workspace.zval)
                    self.zlabel.txt = "Z";
                else
                    self.zlabel.txt = "Ellipsoid Index";
                end
            end

            %%%%
            %%%% Setup the gramian values.
            %%%%

            if  numel(self.workspace.gramian) == 0
                % The gramians of all ellipsoids are equal and identity matrix.
                self.workspace.gramian = repmat(eye(self.ndim, self.ndim), 1, 1, self.nell);
            elseif  numel(self.workspace.gramian) == 1
                % The gramians of all ellipsoids are equal and a multiple of the identity matrix.
                self.workspace.gramian = repmat(self.workspace.gramian * eye(self.ndim, self.ndim), 1, 1, self.nell);
            elseif  size(self.workspace.gramian, 1) == self.ndim && size(self.workspace.gramian, 2) == self.ndim && size(self.workspace.gramian, 3) == 1
                % All ellipsoids have the same user-prescribed self.workspace.gramian matrix.
                self.workspace.gramian = repmat(self.workspace.gramian(:, :, 1), 1, 1, self.nell);
            elseif  size(self.workspace.gramian, 1) == 1 && size(self.workspace.gramian, 2) == 1 && size(self.workspace.gramian, 3) == self.nell
                % Each ellipse has its unique diagonal self.workspace.gramian matrix whose diagonal values are all equal and prescribed by the user.
                temp = zeros(self.ndim, self.ndim, self.nell);
                for iell = 1 : self.nell
                    temp(:, :, iell) = self.workspace.gramian(1, 1, iell) * eye(self.ndim, self.ndim);
                end
                self.workspace.gramian = temp;
            elseif  size(self.workspace.gramian, 1) ~= self.ndim && size(self.workspace.gramian, 2) ~= self.ndim && size(self.workspace.gramian, 3) ~= self.nell
                help("pm.vis.SubplotEllipse3")
                disp("size(gramian)")
                disp( size(self.workspace.gramian) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``gramian`` is incompatible with the" + newline ...
                        + "inferred number of ellipsoid dimensions ``ndim`` and the number of ellipsoids ``nell``." + newline ...
                        + "For more information, see the documentation of the ``gramian`` component displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Setup the center values.
            %%%%

            if  numel(self.workspace.center) == 0
                % The centers of all ellipsoids have the same coordinates and all coordinate axes have zero values.
                self.workspace.center = zeros(self.ndim, self.nell);
            elseif  numel(self.workspace.center) == 1
                % The centers of all ellipsoids have the same coordinates and all coordinate axes have equal non-zero values.
                self.workspace.center = self.workspace.center * ones(self.ndim, self.nell);
            elseif  size(self.workspace.center, 1) == self.ndim && size(self.workspace.center, 2) == 1
                % The centers of all ellipsoids have the same coordinates.
                self.workspace.center = repmat(self.workspace.center, 1, self.nell);
            elseif  size(self.workspace.center, 1) == 1 && size(self.workspace.center, 2) == self.nell
                % All coordinates of each ellipse self.workspace.center are the same.
                self.workspace.center = repmat(self.workspace.center, self.ndim, 1);
            elseif  size(self.workspace.center, 1) ~= self.ndim || size(self.workspace.center, 2) ~= self.nell
                help("pm.vis.SubplotEllipse3")
                disp("size(self.workspace.center)")
                disp( size(self.workspace.center) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``center`` is incompatible with the" + newline ...
                        + "inferred number of ellipsoid dimensions ``ndim`` and the number of ellipsoids ``nell``." + newline ...
                        + "For more information, see the documentation of the ``center`` component displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Setup the zval values.
            %%%%

            if  numel(self.workspace.zval) == 0
                % 2D plot in 3D visualization at zero z value.
                self.workspace.zval = repmat(1 : self.nell, self.npnt, 1);
                if  isempty(self.zscale)
                    self.zscale = "log";
                end
            elseif  numel(self.workspace.zval) == 1
                % 2D plot in 3D visualization at a non-zero z value.
                self.workspace.zval = self.workspace.zval * ones(self.npnt, self.nell);
            elseif  size(self.workspace.zval, 1) == self.npnt && size(self.workspace.zval, 2) == 1
                % The corresponding points on all ellipsoids have the same z values.
                self.workspace.zval = repmat(self.workspace.zval, 1, self.nell);
            elseif  size(self.workspace.zval, 1) == 1 && size(self.workspace.zval, 2) == self.nell
                % All points on each ellipse have the same z values.
                self.workspace.zval = repmat(self.workspace.zval, self.npnt, 1);
            elseif  size(self.workspace.zval, 1) ~= self.npnt || size(self.workspace.zval, 2) ~= self.nell
                help("pm.vis.SubplotEllipse3")
                disp("size(self.workspace.zval)")
                disp( size(self.workspace.zval) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``zval`` is incompatible with" + newline ...
                        + "the inferred number of ellipsoids ``nell`` in the input self.workspace." + newline ...
                        + "For more information, see the documentation of the ``zval`` component displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Setup the cval values. This is relevant only to colormap mode.
            %%%%

            if  numel(self.workspace.cval) == 0
                self.workspace.cval = repmat(1 : self.nell, self.npnt, 1);
                if ~pm.array.len(self.axes.colorScale)
                    self.axes.colorScale = "log";
                end
            elseif  numel(self.workspace.cval) == 1
                % 2D plot in 3D visualization at a non-zero z value.
                self.workspace.cval = self.workspace.cval * ones(self.npnt, self.nell);
            elseif  size(self.workspace.cval, 1) == self.npnt && size(self.workspace.cval, 2) == 1
                % The corresponding points on all ellipsoids have the same z values.
                self.workspace.cval = repmat(self.workspace.cval, 1, self.nell);
            elseif  size(self.workspace.cval, 1) == 1 && size(self.workspace.cval, 2) == self.nell
                % All points on each ellipse have the same z values.
                self.workspace.cval = repmat(self.workspace.cval, self.npnt, 1);
            elseif  size(self.workspace.cval, 1) ~= self.npnt || size(self.workspace.cval, 2) ~= self.nell
                help("pm.vis.SubplotEllipse3")
                disp("size(self.workspace.cval)")
                disp( size(self.workspace.cval) )
                disp("[self.ndim, self.nell]")
                disp( [self.ndim, self.nell] )
                error   ( newline ...
                        + "The shapes of the specified ``cval`` is incompatible with" + newline ...
                        + "the inferred number of ellipsoids ``nell`` in the input self.workspace." + newline ...
                        + "For more information, see the documentation of the ``cval`` component displayed above." + newline ...
                        + newline ...
                        );
            end

           %%%%
           %%%% Check the consistency of shape of ``self.ellindex``.
           %%%%

            if ~isempty(self.ellindex)
                self.workspace.ellindex = self.ellindex(:);
            else
                self.workspace.ellindex = pm.array.logrange(1, self.nell, 75);
            end

            %%%%
            %%%% Generate column indices.
            %%%%

            fields = ["colx", "coly", "colz", "colc"];
            for icol = 1 : length(fields)
                field = fields(icol);
                self.(field) = length(fields) * ((1 : length(self.workspace.ellindex) * numel(self.workspace.dimx)) - 1) + icol;
            end

            %%%%
            %%%% Generate ellipse self.workspace.
            %%%%

            gramtmp = zeros(2, 2);
            bcrd = zeros(size(self.workspace.zval, 1), length(fields) * length(self.workspace.ellindex) * numel(self.workspace.dimx));
            for idim = 1 : numel(self.workspace.dimx)
                for iell = 1 : length(self.workspace.ellindex)
                    jell = self.workspace.ellindex(iell);
                    icol = (iell - 1) * length(fields) + 1;
                    iset = (idim - 1) * length(fields) * length(self.workspace.ellindex);
                    dimxy = [self.workspace.dimx(idim), self.workspace.dimy(idim)];
                    if  self.workspace.dimx(idim) ~= self.workspace.dimy(idim)
                        gramtmp(:, :) = self.workspace.gramian(dimxy, dimxy, jell);
                    else
                        sigmasq = self.workspace.gramian(self.workspace.dimx(idim), self.workspace.dimx(idim), jell);
                        gramtmp(:, :) = [sigmasq, 0; 0, sigmasq];
                    end
                    bcrd(:, iset + icol : iset + icol + 1) = pm.geom.ell2.getBorder ( gramtmp ...
                                                                                    , self.workspace.center(dimxy, jell) ...
                                                                                    , size(self.workspace.zval, 1) ...
                                                                                    );
                    bcrd(:, iset + icol + 2) = self.workspace.zval(:, jell);
                    if  3 < length(fields)
                        bcrd(:, iset + icol + 3) = self.workspace.cval(:, jell);
                    end
                end
            end

            %%%%
            %%%% Generate the dataframe.
            %%%%

            self.df = pm.data.DataFrame(bcrd);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Generate a plot from the selected dimensions
        %>  of the input ellipsoid data to the object constructor.<br>
        %>
        %>  \note
        %>  This method is pure, except for the changes in ``fout`` component.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.SubplotEllipse3](@ref SubplotEllipse3)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      pm.vis.SubplotEllipse3.make();
        %>      pm.vis.SubplotEllipse3.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>  \include{lineno} example/vis/SubplotEllipse3/main.m
        %>  \vis{make}
        %>  <br>\image html example/vis/SubplotEllipse3/SubplotEllipse3.1.png width=700
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 5:46 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            make@pm.vis.SubplotLineScatter3(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if  self.colormap.enabled && self.colorbar.enabled
                if  prod(self.workspace.shape.cval) < 1
                    % input cval is empty
                    %self.fout.colorbar.Ticks = self.fout.axes.ZTick;
                    %%self.fout.colorbar.Limits = self.fout.axes.ZLim;
                    ylabel(self.fout.colorbar, "Ellipsoid Index");
                else
                    ylabel(self.fout.colorbar, "Color Data");
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end