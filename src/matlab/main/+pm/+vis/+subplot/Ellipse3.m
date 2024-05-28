%>  \brief
%>  This is the Ellipse3 class for generating
%>  instances of 3-dimensional Ellipse3 plots
%>  based on the relevant MATLAB
%>  intrinsic functions.
classdef Ellipse3 < pm.vis.subplot.LineScatter3
    properties(Access = public)
        %>
        %>  \param  names   :   The vector of MATLAB strings containing the names of dimensions
        %>                      of the space within which the Gramian matrices are defined.
        %>                      The default is ``"Dimension i"`` where ``i`` is replaced
        %>                      by the ID of the corresponding axis.
        %>
        names = [];
        %>
        %>  \param  dims        :   The vector of size ``2`` of MATLAB whole numbers representing the dimensions
        %>                          of the input ``gramian`` and ``center`` data to use for 2D ellipsoid visualization.
        %>                          The default is ``[1, 2]``.
        %>
        dims = [];
        %>
        %>  \param  ellindex    :   The vector of MATLAB whole numbers representing the ellindex of the 2D
        %>                          ellipsoids to display, represented by the input ``gramian`` and ``center``.
        %>                          The specified ellindex will serve as the visualization values on the z-axis.
        %>                          The default is ``ellindex = pm.array.logrange(start = 1, stop = nell, count = 75);``,
        %>                          where ``nell`` represents the number of ellipsoids to visualize.
        %>
        ellindex = [];
        %>
        %>  \param  gramian     :   A scalar object of class ``pm.data.DataRef``
        %>                          containing the user-specified Gramian data to visualize.
        %>
        gramian = [];
        %>
        %>  \param  center      :   A scalar object of class ``pm.data.DataRef``
        %>                          containing the user-specified center data to visualize.
        %>
        center = [];
        %>
        %>  \param  zval        :   A scalar object of class ``pm.data.DataRef``
        %>                          containing the user-specified z-axis data to visualize.
        %>
        zval = [];
        %>
        %>  \param  cval        :   A scalar object of class ``pm.data.DataRef``
        %>                          containing the user-specified color data in colormap.
        %>                          If empty or unspecified, it will be set to ``zval`` component.
        %>
        cval = [];
    end

    properties(Access = public, Hidden)
        %>
        %>  \param  ndim        :   A scalar MATLAB whole number representing the total number of the
        %>                          dimensions of ellipsoids identified based on the input that can be visualized.
        %>                          It is set to ``max(2, size(gramian, 2), size(center, 1))``.
        %>
        ndim = [];
        %>
        %>  \param  nell        :   A scalar MATLAB whole number representing the total number of
        %>                          ellipsoids identified based on the input that can be visualized.
        %>                          It is set to ``max(20, size(gramian, 3), size(center, 2), size(zval, 2))``.
        %>
        nell = [];
        %>
        %>  \param  npnt        :   The scalar MATLAB whole number representing the number of
        %>                          points with which the 2D ellipsoid boundaries are delineated.
        %>                          It is set to ``max(size(zval, 1), size(cval, 1))`` or otherwise, ``100``.
        %>
        npnt = [];
        %>
        %>  \param  workspace   :   A scalar MATLAB ``struct()`` containing the current visualization information.
        %>                          This information is updated with every call to the ``make()`` method.
        %>                          The contents of this component are temporary and must remain read-only.
        %>                          This component is solely provided for better insight into the internal
        %>                          workings of the ``make()`` method of the parent object and
        %>                          the resulting visualization.
        %>
        workspace = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>
        %>  \details
        %>  In the following documentation,
        %>  <ol>
        %>      <li>    the variable ``ndim`` represents the number of ellipsoids in the input data.
        %>              The value of ``ndim`` is inferred from the shapes
        %>              of the input ``gramian`` and ``center`` arguments
        %>              as ``max([2, size(gramian, 1), size(center, 1)])``.
        %>
        %>      <li>    the variable ``nell`` represents the number of ellipsoids in the input data.
        %>              The value of ``nell`` is inferred from the shapes
        %>              of the input ``gramian``, ``center``, ``zval``, and ``cval`` arguments
        %>              as ``max([size(gramian, 3), size(center, 2), size(zval, 2), size(cval, 2)])``.
        %>              If the above expression yields zero, then ``nell`` is set to ``75``.
        %>
        %>      <li>    the variable ``npnt`` represents the number of points used in visualizing each ellipsoid.
        %>              The value of ``npnt`` is inferred from the shapes of the input arguments
        %>              ``zval`` and ``czval`` as ``max(size(zval, 1), size(cval, 1))``.
        %>              If the above expression yields zero or one, then ``npnt`` is set to ``100``.
        %>  </ol>
        %>
        %>  \note
        %>  To generate a 2D plot of ellipsoids, simply execute the
        %>  MATLAB command ``view(2)`` to change the camera view to 2D.
        %>
        %>  \param[in]  gramian     :   The MATLAB (function handle returning an) object that can be either<br>
        %>                          <ol>
        %>                              <li>    an empty object (e.g., ``[]``), in which case,
        %>                                      the value ``repmat(eye(ndim, ndim), 1, 1, nell)`` will be used.
        %>
        %>                              <li>    a scalar, containing the value of the diagonal elements
        %>                                      of the diagonal matrices representing the Gramian matrices
        %>                                      of all ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a vector of shape ``[1, 1, nell]`` such that ``gramian(1, 1, iell)``
        %>                                      represents the diagonal value to be used in constructing the Gramian
        %>                                      matrix of the ellipsoid ``iell` in the input ``nell`` ellipsoid data.
        %>
        %>                              <li>    a matrix of shape ``[ndim, ndim, 1]`` used in constructing the Gramian
        %>                                      matrices of all ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a tensor of shape ``[ndim, ndim, nell]`` such that ``gramian(:, :, iell)``
        %>                                      represents the Gramian matrix of the ellipsoid ``iell` in the input data.
        %>                          </ol>
        %>
        %>  \note
        %>  While it is possible to pass a data directly to this class,
        %>  it is highly recommended to pass a function handle that returns
        %>  such data when called, allowing the visualization data to be
        %>  dynamically updated when needed.
        %>
        %>  \param[in]  center      :   The MATLAB (function handle returning an) object that can be either
        %>                          <ol>
        %>                              <li>    an empty object (e.g., ``[]``), in which case,
        %>                                      the value ``zeros(ndim, nell)`` will be used.
        %>
        %>                              <li>    a scalar, containing the value to assign to all coordinates
        %>                                      of the centers of all ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a vector of shape ``[1, nell]`` such that ``center(1, iell)``
        %>                                      represents the value of all coordinates of the center of the
        %>                                      ellipsoid ``iell` in the input ``nell`` ellipsoid data.
        %>
        %>                              <li>    a vector of shape ``[ndim, 1]`` such that ``center(idim, 1)``
        %>                                      represents the value of ``idim`` coordinate of the centers of
        %>                                      all ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a matrix of shape ``[ndim, nell]`` such that ``center(:, iell)`` represents
        %>                                      the coordinates of the center of the ellipsoid ``iell` in the input data.
        %>                          </ol>
        %>
        %>  \note
        %>  While it is possible to pass data directly to this class,
        %>  it is highly recommended to pass a function handle that returns
        %>  such data when called, allowing the visualization data to be
        %>  dynamically updated when needed.
        %>
        %>  \param[in]  zval        :   The MATLAB (function handle returning an) object that can be either
        %>                          <ol>
        %>                              <li>    an empty object (e.g., ``[]``), in which case,
        %>                                      the value ``repmat(1 : nell, npnt, 1)`` will be used.
        %>
        %>                              <li>    a scalar, containing the value to assign to the z-axis
        %>                                      coordinates of all ``npnt`` points used in the visualization
        %>                                      of ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a vector of shape ``[1, nell]`` such that ``zval(1, iell)``
        %>                                      represents the value of z-axis coordinates of all ``npnt``
        %>                                      points used in the visualization of the ellipsoid
        %>                                      ``iell` in the input ``nell`` ellipsoid data.
        %>
        %>                              <li>    a vector of shape ``[npnt, 1]`` such that ``zval(ipnt, 1)``
        %>                                      contains the z-axis value of ``ipnt`` points in the representations
        %>                                      of all ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a matrix of shape ``[npnt, nell]`` such that ``zval(:, iell)`` represents
        %>                                      the z-axis coordinates of all ``npnt`` points used in the visualization of
        %>                                      the ellipsoid ``iell` in the input data.
        %>                          </ol>
        %>
        %>  \note
        %>  While it is possible to pass data directly to this class,
        %>  it is highly recommended to pass a function handle that returns
        %>  such data when called, allowing the visualization data to be
        %>  dynamically updated when needed.
        %>
        %>  \param[in]  cval        :   The MATLAB (function handle returning an) object that can be either
        %>                          <ol>
        %>                              <li>    an empty object (e.g., ``[]``), in which case,
        %>                                      the input value ``zval`` will be used.
        %>
        %>                              <li>    a scalar, containing the value to assign to color values
        %>                                      of all ``npnt`` points used in the visualization
        %>                                      of  ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a vector of shape ``[1, nell]`` such that ``cval(1, iell)``
        %>                                      represents the color values of all ``npnt`` points used in the
        %>                                      visualization of the ellipsoid ``iell` in the input ``nell`` ellipsoid data.
        %>
        %>                              <li>    a vector of shape ``[npnt, 1]`` such that ``cval(ipnt, 1)``
        %>                                      contains the color values of ``ipnt`` points in the representations
        %>                                      of all ``nell`` ellipsoids in the input data.
        %>
        %>                              <li>    a matrix of shape ``[npnt, nell]`` such that ``center(:, iell)`` represents
        %>                                      the color values of all ``npnt`` points used in the visualization of
        %>                                      the ellipsoid ``iell` in the input data.
        %>                          </ol>
        %>
        %>  \note
        %>  While it is possible to pass data directly to this class,
        %>  it is highly recommended to pass a function handle that returns
        %>  such data when called, allowing the visualization data to be
        %>          dynamically updated when needed.
        %>
        %>  \note
        %>  The input value for ``cval`` is used only if ``colormap.enabled``
        %>  component of the output object of class ``pm.vis.subplot.Ellipse3``
        %>  is set to ``true``.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  See below and also the documentation of the attributes
        %>  of the superclass ``pm.vis.subplot.LineScatter3``.
        %>
        %>  \return
        %>  An object of class ``pm.vis.subplot.Ellipse3``.
        %>
        %>  \interface{Ellipse3}
        %>  \code{.m}
        %>      s = pm.vis.subplot.Ellipse3();
        %>      s = pm.vis.subplot.Ellipse3(gramian);
        %>      s = pm.vis.subplot.Ellipse3(gramian, center);
        %>      s = pm.vis.subplot.Ellipse3(gramian, center, zval);
        %>      s = pm.vis.subplot.Ellipse3(gramian, center, zval, cval);
        %>      s = pm.vis.subplot.Ellipse3(gramian, center, zval, cval, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Ellipse3}
        %>
        %>      s = pm.vis.subplot.Ellipse3();
        %>      s.make("dims", [1, 2]);
        %>
        %>  \final{Ellipse3}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 5:35 PM, University of Texas at Arlington<br>
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

            varargin =  { "gramian", pm.data.DataRef(gramian) ...
                        , "center", pm.data.DataRef(center) ...
                        , "zval", pm.data.DataRef(zval) ...
                        , "cval", pm.data.DataRef(cval) ...
                        , varargin{:} ...
                        };
            self = self@pm.vis.subplot.LineScatter3([], varargin{:})

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>  \brief
        %>  Reset the properties of the plot to the original default settings.
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.
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
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.subplot.Ellipse3.reset() # reset the plot to the default settings.
        %>
        %>  \endcode
        %>  \final{len}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 5:42 PM, University of Texas at Arlington<br>
        %>
        function reset(self, varargin)
            self.dims = [];
            self.npnt = [];
            self.names = [];
            self.ellindex = [];
            reset@pm.vis.subplot.LineScatter3(self, varargin{:});
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
        %>  Prepare the subplot for visualization.
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``premake()`` method.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      pm.vis.subplot.Ellipse3.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{premake}
        %>
        %>      s = pm.vis.subplot.Ellipse3();
        %>      s.premake("dims", [1, 2], "colormap", {"map", "autumn"})
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 5:45 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
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

            premake@pm.vis.subplot.LineScatter3(self);

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
                help("pm.vis.subplot.Ellipse3")
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
            %%%% Set the default dims or verify the custom values to visualize.
            %%%%

            if  isempty(self.dims)
                self.dims = [1, 2];
            elseif  size(self.dims, 2) ~= 2
                help("pm.vis.subplot.Ellipse3")
                disp("self.dims")
                disp( self.dims )
                disp("size(self.dims)")
                disp( size(self.dims) )
                error   ( newline ...
                        + "The size of the specified ``dims`` must be ``2``." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end

           %%%%
           %%%% Check the consistency of shape of ``self.names``.
           %%%%

            if ~isempty(self.names) && pm.array.len(self.names) ~= self.ndim
                help("pm.vis.subplot.Ellipse3")
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
                self.xlabel.txt = string(self.workspace.names(self.dims(:, 1)));
            end
            if ~pm.array.len(self.ylabel.txt)
                self.ylabel.txt = string(self.workspace.names(self.dims(:, 2)));
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
                help("pm.vis.subplot.Ellipse3")
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
                help("pm.vis.subplot.Ellipse3")
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
                help("pm.vis.subplot.Ellipse3")
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
                help("pm.vis.subplot.Ellipse3")
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
                self.(field) = length(fields) * ((1 : length(self.workspace.ellindex) * size(self.dims, 1)) - 1) + icol;
            end

            %%%%
            %%%% Generate ellipse self.workspace.
            %%%%

            bcrd = zeros(size(self.workspace.zval, 1), length(fields) * length(self.workspace.ellindex) * size(self.dims, 1));
            for idim = 1 : size(self.dims, 1)
                for iell = 1 : length(self.workspace.ellindex)
                    jell = self.workspace.ellindex(iell);
                    icol = (iell - 1) * length(fields) + 1;
                    iset = (idim - 1) * length(fields) * length(self.workspace.ellindex);
                    bcrd(:, iset + icol : iset + icol + 1) = pm.geom.ell2.getBorder ( self.workspace.gramian(self.dims(idim, :), self.dims(idim, :), jell) ...
                                                                                    , self.workspace.center(self.dims(idim, :), jell) ...
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
        %>
        %>      pm.vis.subplot.Ellipse3.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>
        %>      s = pm.vis.subplot.Ellipse3();
        %>      s.make("dims", [1, 2]);
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 5:46 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            make@pm.vis.subplot.LineScatter3(self, varargin{:});

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