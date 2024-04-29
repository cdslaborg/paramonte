classdef Axes < pm.matlab.Handle
    %
    %   This is the abstract class for generating instances of objects
    %   that contain the specifications of various types of plots.
    %
    %   This class primarily serves as the superclass for
    %   the visualization-ready subclass ``pm.vis.subplot.Subplot``
    %   and its subclasses, all accessible to the end users.
    %
    %   Parameters
    %   ----------
    %
    %       ptype
    %
    %           The input scalar MATLAB string containing the name of the
    %           subclass that whose parent is Axes (e.g., "heatmap").
    %           Supported plot names are:
    %
    %               line
    %               line3
    %               scatter
    %               scatter3
    %               lineScatter
    %               lineScatter3
    %               histogram2
    %               histogram
    %               contour3
    %               contourf
    %               contour
    %               histfit
    %               heatmap
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``premake()`` method.
    %
    %   Returns
    %   -------
    %
    %       self
    %
    %           The output scalar object of class ``pm.vis.Axes``.
    %
    %   Interface
    %   ---------
    %
    %       axes = pm.vis.Axes(ptype);
    %
    %   Attributes
    %   ----------
    %
    %       The following is the list of all class attributes that are dynamically
    %       added to the instantiated class objects based on the specified input plot type.
    %       See also the explicit class and superclass attributes not listed below.
    %
    %       \devnote
    %
    %           While dynamic addition of class attributes is not ideal, the current
    %           design was deemed unavoidable and best, given the constraints of the
    %           MATLAB language and visualization tools.
    %
    %       axes (available for all plots except heatmap)
    %
    %           A MATLAB ``struct`` whose fields and values are passed as
    %           keyword arguments to the MATLAB intrinsic ``set()`` for
    %           the current active axes object in the plot ``gca()``.
    %
    %       colorbar (available for all axes types that allow color-mapping)
    %
    %           A MATLAB ``struct`` whose fields and their values will
    %           be passed as keyword arguments to the MATLAB intrinsic ``colorbar``.
    %           The following are the default components of ``colorbar``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   colorbar will be applied to the axes.
    %
    %               others
    %
    %                   See the acceptable keyword arguments of the MATLAB intrinsic ``colorbar()``.
    %
    %           Example usage:
    %
    %               self.colorbar.enabled = true;
    %               self.colorbar.location = "west";
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %               For example, ``colorbar.color`` and ``colorbar.Color`` are the same,
    %               and only one of the two will be processed.
    %
    %       colormap (available for all axes types that allow color-mapping)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``colormap``.
    %           The following are the default components of ``colormap``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   colormap will be applied to the axes.
    %
    %               map
    %
    %                   A string or a vector of color triplets or any other value
    %                   that the intrinsic MATLAB ``colormap`` accepts as input.
    %
    %           This option is relevant only to visualizations that allow color-mapping.
    %
    %           Example usage:
    %
    %               1.  self.colormap.enabled = true;
    %               2.  self.colormap.map = "winter";
    %               2.  self.colormap.map = "winter";
    %               3.  self.colormap.map = 'default';
    %               4.  self.colormap.map = pm.vis.cmap.redblue();
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %               For example, ``colormap.map`` and ``colormap.Map`` are the same,
    %               and only one of the two will be processed.
    %
    %       contour (available only for contour axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``contour``.
    %           The following are the default components of ``contour``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   contour will be added to the axes.
    %
    %               levels
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.
    %
    %               lineSpec
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.
    %
    %               others
    %
    %                   See the acceptable keyword arguments of the MATLAB intrinsic ``contour()``.
    %
    %           Example usage:
    %
    %               1.  self.contour.enabled = true;
    %               2.  self.contour.lineWidth = "none";
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       contour3 (available only for contour3 axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``contour3``.
    %           The following are the default components of ``contour3``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   contour3 will be added to the axes.
    %
    %               levels
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.
    %
    %               lineSpec
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.
    %
    %               others
    %
    %                   See the acceptable keyword arguments of the MATLAB intrinsic ``contour3()``.
    %
    %           Example usage:
    %
    %               1.  self.contour3.enabled = true;
    %               2.  self.contour3.lineWidth = "none";
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       contourf (available only for contourf axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``contourf``.
    %           The following are the default components of ``contourf``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   contourf will be added to the axes.
    %
    %               levels
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.
    %
    %               lineSpec
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``contourf()``.
    %
    %               others
    %
    %                   See the acceptable keyword arguments of the MATLAB intrinsic ``contourf()``.
    %
    %           Example usage:
    %
    %               1.  self.contourf.enabled = true;
    %               2.  self.contourf.lineWidth = "none";
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       histfit (available only for histfit axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``histfit``.
    %           The following are the default components of ``histfit``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   histfit will be added to the axes.
    %
    %               nbins
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histfit()``.
    %
    %               dist
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histfit()``.
    %
    %           Example usage:
    %
    %               1.  self.histfit.enabled = true;
    %               2.  self.histfit.nbins = 20;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       histogram (available only for histogram axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``histogram``.
    %           The following are the default components of ``histogram``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   histogram will be added to the axes.
    %
    %               nbins
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histogram()``.
    %
    %               edges
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histogram()``.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``histogram()``.
    %
    %           Example usage:
    %
    %               1.  self.histogram.enabled = true;
    %               2.  self.histogram.edgeColor = "none";
    %               3.  self.histogram.nbins = 20;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       histogram2 (available only for histogram2 axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``histogram2``.
    %           The following are the default components of ``histogram2``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   histogram2 will be added to the axes.
    %
    %               nbins
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histogram2()``.
    %
    %               xedges
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histogram2()``.
    %
    %               yedges
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``histogram2()``.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``histogram2()``.
    %
    %           Example usage:
    %
    %               1.  self.histogram2.enabled = true;
    %               2.  self.histogram2.edgeColor = "none";
    %               3.  self.histogram2.nbins = 20;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       legend (available for all axes types except heatmap)
    %
    %           A MATLAB ``struct`` whose fields and values are passed
    %           as keyword arguments to the MATLAB intrinsic ``title()``.
    %
    %       maxnoise (available only for contour/contourf/contour3 axes types)
    %
    %           A float indicating the threshold below which the kernel density
    %           estimate is considered to be noise and is rounded to zero.
    %           The higher this value is, the less noise will be
    %           visible in the resulting contour plots.
    %           If empty, the default value is ``0.001``.
    %
    %       plot (available only for line/lineScatter axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``plot``.
    %           The following are the default components of ``plot``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   plot will be added to the axes.
    %
    %               lineSpec
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``plot()``.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``plot()``.
    %
    %           Example usage:
    %
    %               1.  self.plot.enabled = true;
    %               2.  self.plot.lineWidth = 1;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       plot3 (available only for line3/lineScatter3 axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``plot3``.
    %           The following are the default components of ``plot3``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   plot3 will be added to the axes.
    %
    %               lineSpec
    %
    %                   See the corresponding positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``plot3()``.
    %
    %           Example usage:
    %
    %               1.  self.plot3.enabled = true;
    %               2.  self.plot3.lineWidth = 1;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       precision (available only for heatmap axes types)
    %
    %           A scalar integer representing the number of digits after
    %           the decimal point for the values that appear in each cell
    %           of the heatmap. The default value is set by MATLAB.
    %
    %       resolution (available only for contour/contourf/contour3 axes types)
    %
    %           A scalar integer indicating the grid resolution for discretization of
    %           the data during the kernel density estimation. It must be a power of
    %           two, otherwise it will be changed to the next power of two at the
    %           time of using it. If empty, the default value is ``2^9``.
    %
    %       scatter (available only for scatter/lineScatter axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``scatter``.
    %           The following are the default components of ``scatter``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   scatter will be added to the axes.
    %
    %               size
    %
    %                   See the corresponding ``sz`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               color
    %
    %                   See the corresponding ``C`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               filled
    %
    %                   See the corresponding ``filled`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               marker
    %
    %                   See the corresponding ``mkr`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``scatter()``.
    %
    %           Example usage:
    %
    %               1.  self.scatter.enabled = true; % add scatter()
    %               4.  self.scatter.color = "red"; % set the points color
    %               2.  self.scatter.marker = "."; % set the marker type
    %               3.  self.scatter.size = 10; % set the point size
    %               4.  self.scatter.lineWidth = 0;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       scatter3 (available only for scatter3/lineScatter3 axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``scatter3``.
    %           The following are the default components of ``scatter3``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   scatter3 will be added to the axes.
    %
    %               size
    %
    %                   See the corresponding ``sz`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               color
    %
    %                   See the corresponding ``C`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               filled
    %
    %                   See the corresponding ``filled`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               marker
    %
    %                   See the corresponding ``mkr`` positional argument of the MATLAB intrinsic ``plot3()``.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``scatter3()``.
    %
    %           Example usage:
    %
    %               1.  self.scatter3.enabled = true; % add scatter3()
    %               4.  self.scatter3.color = "red"; % set the points color
    %               2.  self.scatter3.marker = "."; % set the marker type
    %               3.  self.scatter3.size = 10; % set the point size
    %               4.  self.scatter3.lineWidth = 0;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       surface (available only for line/lineScatter/line3/lineScatter3 axes types)
    %
    %           A MATLAB ``struct`` whose fields and their values will be passed
    %           as keyword arguments to the MATLAB intrinsic ``surface``.
    %           The following are the default components of ``surface``:
    %
    %               enabled
    %
    %                   A logical value. If ``true``, the
    %                   surface will be added to the axes.
    %
    %               others
    %
    %                   See the corresponding acceptable keyword arguments of the MATLAB intrinsic ``surface()``.
    %
    %           Example usage:
    %
    %               1.  self.surface.enabled = true;
    %               2.  self.surface.edgeColor = 1;
    %
    %           \warning
    %
    %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
    %               Hence, ensure you do not add the same keyword as multiple different fields.
    %
    %       target (available only for all axes types except heatmap)
    %
    %           An object of class ``pm.vis.Target`` for adding target values to the plots.
    %           For more information, see the documentation for the class ``pm.vis.Target``.
    %
    %       title (available for all axes types)
    %
    %           A MATLAB ``struct`` whose fields and values are passed
    %           as keyword arguments to the MATLAB intrinsic ``title()``.
    %
    %       xlabel (available for all axes types)
    %
    %           A MATLAB ``struct`` whose fields and values are passed
    %           as keyword arguments to the MATLAB intrinsic ``xlabel()``.
    %
    %       ylabel (available for all axes types)
    %
    %           A MATLAB ``struct`` whose fields and values are passed
    %           as keyword arguments to the MATLAB intrinsic ``ylabel()``.
    %
    %       zlabel (available only for all tri-axes axes types)
    %
    %           A MATLAB ``struct`` whose fields and values are passed
    %           as keyword arguments to the MATLAB intrinsic ``zlabel()``.
    %
    %       xlim (available for all axes types)
    %
    %           A MATLAB vector of length ``2`` whose fields and values are
    %           passed as keyword arguments to the MATLAB intrinsic ``xlim()``.
    %
    %       ylim (available for all axes types)
    %
    %           A MATLAB vector of length ``2`` whose fields and values are
    %           passed as keyword arguments to the MATLAB intrinsic ``ylim()``.
    %
    %       zlim (available only for all tri-axes axes types)
    %
    %           A MATLAB vector of length ``2`` whose fields and values are
    %           passed as keyword arguments to the MATLAB intrinsic ``zlim()``.
    %
    %       xscale (available for all axes types)
    %
    %           A MATLAB string whose value is passed directly to the MATLAB intrinsic
    %           ``xscale()`` to set the axis scale to either logarithmic or linear.
    %           Possible values are: ``"log"``, ``"linear"``.
    %           The default behavior is set by MATLAB.
    %
    %       yscale (available for all axes types)
    %
    %           A MATLAB string whose value is passed directly to the MATLAB intrinsic
    %           ``yscale()`` to set the axis scale to either logarithmic or linear.
    %           Possible values are: ``"log"``, ``"linear"``.
    %           The default behavior is set by MATLAB.
    %
    %       zscale (available only for all tri-axes axes types)
    %
    %           A MATLAB string whose value is passed directly to the MATLAB intrinsic
    %           ``zscale()`` to set the axis scale to either logarithmic or linear.
    %           Possible values are: ``"log"``, ``"linear"``.
    %           The default behavior is set by MATLAB.
    %
    properties(Access = protected, Hidden)
        type = struct();
        cenabled = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Axes(ptype, varargin)

            if  nargin < 1 || ~pm.introspection.istype(ptype, "string", 1) || ~pm.array.len(ptype)
                help("pm.vis.Axes");
                error   ( newline ...
                        + "The input argument ``ptype`` is missing." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            else
                % lower the first character.
                ptype = string(ptype);
                ptype{1} = lower(ptype{1});
            end

            self.type = struct();
            self.type.name = ptype;
            pnamel = lower(ptype);
            self.type.is.line           = strcmpi(pnamel, "line"        );
            self.type.is.line3          = strcmpi(pnamel, "line3"       );
            self.type.is.scatter        = strcmpi(pnamel, "scatter"     );
            self.type.is.scatter3       = strcmpi(pnamel, "scatter3"    );
            self.type.is.lineScatter    = strcmpi(pnamel, "lineScatter" );
            self.type.is.lineScatter3   = strcmpi(pnamel, "lineScatter3");
            self.type.is.histogram2     = strcmpi(pnamel, "histogram2"  );
            self.type.is.histogram      = strcmpi(pnamel, "histogram"   );
            self.type.is.contour3       = strcmpi(pnamel, "contour3"    );
            self.type.is.contourf       = strcmpi(pnamel, "contourf"    );
            self.type.is.contour        = strcmpi(pnamel, "contour"     );
            self.type.is.heatmap        = strcmpi(pnamel, "heatmap"     );
            self.type.is.histfit        = strcmpi(pnamel, "histfit"     );
            self.type.is.diffusion      = self.type.is.contour || self.type.is.contourf || self.type.is.contour3;
            self.type.is.density        = self.type.is.diffusion || self.type.is.histfit || self.type.is.histogram || self.type.is.histogram2;
            self.type.is.d3             = self.type.is.line3 || self.type.is.scatter3 || self.type.is.lineScatter3;
            self.type.is.d1             = self.type.is.histfit || self.type.is.histogram || self.type.is.heatmap;
            self.type.is.d2             = ~(self.type.is.d1 || self.type.is.d3);
            self.type.is.triaxes        = self.type.is.d3 || self.type.is.histogram2 || self.type.is.contour3;
            self.type.is.targetable     = false; %xx ~self.type.is.heatmap && ~ self.type.is.d3;
            self.type.has.line          = self.type.is.line || self.type.is.line3 || self.type.is.lineScatter || self.type.is.lineScatter3;
            self.type.has.scatter       = self.type.is.scatter || self.type.is.scatter3 || self.type.is.lineScatter || self.type.is.lineScatter3;

            self.reset(varargin{:}); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self, varargin)
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot and
            %   you want to clean up and go back to the default settings.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.Axes.reset() # reset the plot to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%% axes

            if ~self.type.is.heatmap
                self.newprop("axes", struct());
                self.axes.box = [];
                self.axes.color = [];
                self.axes.colorScale = [];
                self.axes.fontName = [];
                self.axes.fontSize = [];
                self.axes.fontSizeMode = [];
                self.axes.fontSmoothing = [];
                self.axes.fontWeight = [];
                self.axes.xgrid = [];
                self.axes.ygrid = [];
                if  self.type.is.triaxes
                    self.axes.zgrid = [];
                end
                self.axes.enabled = true;
            end

            %%%% title

            self.newprop("title", struct());
            if ~self.type.is.heatmap
                self.title.fontSize = [];
                self.title.interpreter = [];
                self.title.fontWeight = [];
                self.title.color = []; %[0, 0, 0];
            end
            self.title.enabled = [];
            self.title.titletext = [];
            if ~self.type.is.heatmap
                self.title.subtitletext = [];
            end

            %%%% xlabel, xlim

            self.newprop("xlabel", struct());
            if ~self.type.is.heatmap
                self.xlabel.color = []; %[0.15, 0.15, 0.15];
                self.xlabel.fontAngle = [];
                self.xlabel.fontSize = [];
                self.xlabel.fontWeight = [];
                self.xlabel.interpreter = [];
                self.xlabel.rotation = [];
            end
            self.xlabel.enabled = true;
            self.xlabel.txt = [];

            self.newprop("xlim", []);
            self.newprop("xscale", []);

            %%%% ylabel, ylim

            self.newprop("ylim", []);
            self.newprop("ylabel", self.xlabel);
            self.newprop("yscale", []);

            %%%% zlabel, zlim

            if  self.type.is.triaxes
                self.newprop("zlabel", self.xlabel);
                self.newprop("zlim", []);
                self.newprop("zscale", []);
            end

            %%%% colc, colorbar, colormap

            if  self.type.is.heatmap || ~self.type.is.d1

                if  self.type.has.line || self.type.has.scatter
                    self.newprop("colc", {});
                end

                if ~self.type.is.heatmap
                    self.newprop("colorbar", struct());
                    self.colorbar.fontSize = [];
                    self.colorbar.direction = 'normal';
                    self.colorbar.limits = [];
                    self.colorbar.location = 'eastoutside';
                    self.colorbar.ticks = [];
                    self.colorbar.tickLabels = [];
                    self.colorbar.tickLabelInterpreter = 'tex';
                    self.colorbar.enabled = true;
                end

                self.newprop("colormap", struct());
                self.colormap.enabled = true;
                self.colormap.map = []; % 'default'

            end

            %%%% legend

            if ~self.type.is.heatmap
                self.newprop("legend", struct());
                self.legend.box = 'off';
                self.legend.color = 'none';
                self.legend.fontSize = [];
                self.legend.interpreter = 'none';
                self.legend.location = [];
                self.legend.numColumns = [];
                self.legend.textColor = [];
                self.legend.enabled = false;
                self.legend.labels = {};
            end

            %%%% target

            if ~self.type.is.targetable
                self.newprop("target");
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            %%%% heatmap attributes
            %%%%%%%%%%%%%%%%%%%%%%%

            if  self.type.is.heatmap
                self.newprop("heatmap", struct());
                self.newprop("precision", []);
                self.heatmap.colorbarVisible = [];
                self.heatmap.colorLimits = [];
                self.heatmap.missingDataColor = [];
                self.heatmap.fontName = [];
                self.heatmap.fontSize = [];
                self.heatmap.enabled = true;
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            %%%% density attributes
            %%%%%%%%%%%%%%%%%%%%%%%

            if  self.type.is.histfit
                self.newprop("histfit", struct());
                self.histfit.enabled = true;
                self.histfit.dist = [];
                self.histfit.nbins = [];
            end

            if  self.type.is.histogram
                self.newprop("histogram", struct());
                self.histogram.binLimitsMode = [];
                self.histogram.binMethod = [];
                self.histogram.binWidth = [];
                self.histogram.displayStyle = [];
                self.histogram.edgeAlpha = [];
                self.histogram.edgeColor = [];
                self.histogram.faceAlpha = [];
                self.histogram.faceColor = [];
                self.histogram.lineStyle = [];
                self.histogram.lineWidth = [];
                self.histogram.normalization = [];
                self.histogram.orientation = [];
                self.histogram.enabled = true;
                self.histogram.edges = [];
                self.histogram.nbins = [];
            end

            if  self.type.is.histogram2
                self.newprop("histogram2", struct());
                self.histogram2.binWidth = [];
                self.histogram2.xbinLimits = [];
                self.histogram2.xbinLimitsMode = [];
                self.histogram2.ybinLimits = [];
                self.histogram2.ybinLimitsMode = [];
                self.histogram2.binMethod = [];
                self.histogram2.showEmptyBins = [];
                self.histogram2.normalization = [];
                self.histogram2.displayStyle = [];
                self.histogram2.edgeAlpha = [];
                self.histogram2.edgeColor = [];
                self.histogram2.faceAlpha = [];
                self.histogram2.faceColor = [];
                self.histogram2.faceLighting = [];
                self.histogram2.lineStyle = [];
                self.histogram2.lineWidth = [];
                self.histogram2.enabled = true;
                self.histogram2.nbins = [];
                self.histogram2.xedges = [];
                self.histogram2.yedges = [];
            end

            if self.type.is.contour || self.type.is.contourf || self.type.is.contour3
                self.newprop("resolution", 2^9);
                self.newprop("maxnoise", 0.001);
            end

            if  self.type.is.contour
                self.newprop("contour", struct());
                self.contour.color = [];
                self.contour.labelSpacing = [];
                self.contour.lineWidth = [];
                self.contour.showText = [];
                self.contour.enabled = true;
                self.contour.levels = [];
                self.contour.lineSpec = [];
            end

            if  self.type.is.contourf
                self.newprop("contourf", struct());
                self.contourf.edgecolor = [];
                self.contourf.color = [];
                self.contourf.labelSpacing = [];
                self.contourf.lineWidth = [];
                self.contourf.showText = [];
                self.contourf.enabled = true;
                self.contourf.levels = [];
                self.contourf.lineSpec = [];
            end

            if  self.type.is.contour3
                self.newprop("contour3", struct());
                self.contour3.color = [];
                self.contour3.labelSpacing = [];
                self.contour3.lineWidth = [];
                self.contour3.showText = [];
                self.contour3.enabled = true;
                self.contour3.levels = [];
                self.contour3.lineSpec = [];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% line/scatter attributes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if  self.type.is.line || self.type.is.lineScatter
                self.newprop("surface", struct());
                self.newprop("plot", struct());
                self.plot.color = [];
                self.plot.lineWidth = [];
                self.plot.lineStyle = [];
                self.plot.markeredgeColor = [];
                self.plot.markerFaceColor = [];
                self.plot.markerIndices = [];
                self.plot.markerSize = [];
                self.plot.enabled = true;
            end

            if  self.type.is.line3 || self.type.is.lineScatter3
                self.newprop("surface", struct());
                self.newprop("plot3", struct());
                self.plot3.color = {};
                self.plot3.lineWidth = {};
                self.plot3.lineStyle = '-';
                self.plot3.markeredgeColor = [];
                self.plot3.markerFaceColor = [];
                self.plot3.markerIndices = [];
                self.plot3.markerSize = [];
                self.plot3.enabled = true;
            end

            if  isprop(self, "surface")
                self.surface.alignVertexCenters = {};
                self.surface.backFaceLighting = {};
                self.surface.edgeAlpha = {};
                self.surface.edgeColor = {};
                self.surface.edgeLighting = {};
                self.surface.faceAlpha = {};
                self.surface.faceColor = {};
                self.surface.faceLighting = {};
                self.surface.lineStyle = {};
                self.surface.lineWidth = {};
                self.surface.marker = {};
                self.surface.markerEdgeColor = {};
                self.surface.markerFaceColor = {};
                self.surface.markerSize = {};
                self.surface.meshStyle = {};
                self.surface.enabled = [];
            end

            if  self.type.is.scatter || self.type.is.lineScatter
                self.newprop("scatter", struct());
                self.scatter.colorVariable = [];
                self.scatter.lineWidth = [];
                self.scatter.markeredgeColor = [];
                self.scatter.markerFaceColor = [];
                self.scatter.enabled = true;
                self.scatter.color = [];
                self.scatter.filled = [];
                self.scatter.marker = [];
                self.scatter.size = [];
            end

            if  self.type.is.scatter3 || self.type.is.lineScatter3
                self.newprop("scatter3", struct());
                self.scatter3.colorVariable = [];
                self.scatter3.lineWidth = [];
                self.scatter3.markeredgeColor = [];
                self.scatter3.markerFaceColor = [];
                self.scatter3.enabled = true;
                self.scatter3.color = [];
                self.scatter3.filled = [];
                self.scatter3.marker = [];
                self.scatter3.size = [];
            end

            %%%% ensure line plots are monochromatic when axes contains scatter plot.

            if  self.type.has.line
                if  self.type.has.scatter
                    self.surface.enabled = false;
                    if  self.type.is.d2
                        self.plot.color = [];
                    else
                        self.plot3.color = [];
                    end
                else
                    self.surface.enabled = true;
                end
            end

            if  self.type.is.targetable
                self.target = pm.vis.Target();
            end

            self.premake(varargin{:}); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function premake(self, varargin)
            %
            %   Configure the plot settings and specifications and return nothing.
            %
            %   \warning
            %
            %       This method has side-effects by manipulating
            %       the existing attributes of the parent object.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes, before calling the ``premake()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       a = pm.vis.Axes.premake(varargin);
            %
            %   Example
            %   -------
            %
            %       a = pm.vis.Axes(ptype);
            %       a.premake("xlim", [0, 1])
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if ~isempty(varargin)
                self.hash2comp(varargin); % parse arguments
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%% Set heatmap settings.

            if  self.type.is.heatmap
                if ~isempty(self.precision) && ~isnumeric(self.precision)
                    help("pm.vis.Axes");
                    disp("self.precision");
                    disp( self.precision );
                    error   ( newline ...
                            + "The specified precision must be a positive whole number." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end
                self.setKeyVal(self.type.name, "colorbarVisible", "on");
                self.setKeyVal(self.type.name, "missingDataColor", "[0.1500 0.1500 0.1500]");
            end

            %%%% Set histfit/histogram/histogram2/contour/contour3/contourf settings.

            if  self.type.is.histogram
                self.setKeyVal(self.type.name, "edgeColor", "none");
            end

            if  self.type.is.histfit
                self.setKeyVal(self.type.name, "dist", "Normal");
                self.setKeyVal(self.type.name, "nbins", []);
            end

            if  self.type.is.histogram2
                self.setKeyVal(self.type.name, "edgeColor", "none");
                self.setKeyVal(self.type.name, "displayStyle", "bar3");
                self.setKeyVal(self.type.name, "showEmptyBins", "off");
               %self.setKeyVal(self.type.name, "numbins", [100 100]);
                if ~self.colormap.enabled
                    if ~pm.introspection.istype(self.histogram2.faceColor, "string", 1)
                        self.setKeyVal(self.type.name, "faceColor", "auto");
                    elseif self.histogram2.faceColor == "flat"
                        % enforce monochrome by removing the colormapping.
                        self.histogram2.faceColor = "auto";
                    end
                else
                    self.histogram2.faceColor = "flat";
                   %self.setKeyVal(self.type.name, "faceColor", "flat");
                end
            end

            if  self.type.is.diffusion
                if  self.type.is.contour3
                    self.setKeyVal(self.type.name, "color", []);
                    self.setKeyVal(self.type.name, "levels", 50);
                else
                    if  self.type.is.contourf
                        self.setKeyVal(self.type.name, "color", "none");
                        self.setKeyVal(self.type.name, "edgecolor", "none");
                        self.setKeyVal(self.type.name, "levels", 50); % 8
                        %if ~self.colormap.enabled && self.contourf.enabled
                        %    self.colormap.map = flipud(gray);
                        %end
                    else
                        self.setKeyVal(self.type.name, "levels", 10); % 8
                        self.setKeyVal(self.type.name, "color", []);
                    end
                end
                self.setKeyVal("resolution", 2^9);
                self.setKeyVal("maxnoise", 0.001);
                self.setKeyVal(self.type.name, "labelSpacing", 144);
                self.setKeyVal(self.type.name, "showText", "off");
                self.setKeyVal(self.type.name, "lineStyle", "-");
                self.setKeyVal(self.type.name, "lineWidth", 0.5);
            end

            %%%% Set line/scatter settings.

            if  self.type.is.scatter || self.type.is.lineScatter
                self.setKeyVal("scatter", "size", 5);
                self.setKeyVal("scatter", "marker", "o");
                self.setKeyVal("scatter", "filled", true);
                if  self.type.is.lineScatter && (self.scatter.enabled || self.colormap.enabled)
                    self.setKeyVal("plot", "color", uint8([200 200 200 150]));
                end
            end

            if  self.type.is.scatter3 || self.type.is.lineScatter3
                self.setKeyVal("scatter3", "size", 5);
                self.setKeyVal("scatter3", "marker", "o");
                self.setKeyVal("scatter3", "filled", true);
                if  self.type.is.lineScatter3 && (self.scatter3.enabled || self.colormap.enabled)
                    self.setKeyVal("plot3", "color", uint8([200 200 200 150]));
                end
            end

            if  self.type.has.line
                if  self.surface.enabled
                    if self.type.is.d2; self.plot.enabled = false; end
                    if self.type.is.d3; self.plot3.enabled = false; end
                end
                if self.type.is.d2; self.setKeyVal("plot", "lineWidth", 1); end
                if self.type.is.d3; self.setKeyVal("plot3", "lineWidth", 1); end
                self.setKeyVal("surface", "faceColor", "none");
                self.setKeyVal("surface", "edgeColor", "flat");
                self.setKeyVal("surface", "edgeAlpha", 0.5);
                self.setKeyVal("surface", "lineStyle", "-");
                self.setKeyVal("surface", "lineWidth", 1);
                self.setKeyVal("surface", "marker", "none");
            end

            %%%% Set the coloring settings.

            self.cenabled = false;
            self.cenabled = self.cenabled || (self.type.has.scatter && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.has.line && self.surface.enabled && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.is.histogram2 && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.is.diffusion && self.colormap.enabled);
            self.cenabled = self.cenabled || (self.type.is.heatmap && self.colormap.enabled);

            if ~self.cenabled && self.type.has.line
                if self.type.is.d2; self.plot.enabled = true; end
                if self.type.is.d3; self.plot3.enabled = true; end
            end

            if isprop(self, "axes")
                if isfield(self.axes, "box") && isempty(self.axes.box); self.axes.box = "on"; end
                if isfield(self.axes, "xgrid") && isempty(self.axes.xgrid); self.axes.xgrid = "on"; end
                if isfield(self.axes, "ygrid") && isempty(self.axes.ygrid); self.axes.ygrid = "on"; end
                if isfield(self.axes, "zgrid") && isempty(self.axes.zgrid); self.axes.zgrid = "on"; end
            end

            %if  (self.type.is.diffusion || self.type.is.histogram2) && isempty(self.zlabel.txt)
            %    self.zlabel.txt = "Density";
            %end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function hash = comp2hash(self, comp)
            %
            %   Convert the components of the input component ``comp``
            %   of the parent object into a cell array of key-val pairs.
            %
            excludes = {"enabled"};
            if strcmp(comp, "axes")
                excludes = [excludes(:); "labels"; "parent"; "ncol"; "nrow"];
            elseif strcmp(comp, "colorbar")
                excludes = [excludes(:); "width"; "height"];
            elseif strcmp(comp, "colormap")
                excludes = [excludes(:); "map"];
            elseif strcmp(comp, "contour3") || strcmp(comp, "contourf") || strcmp(comp, "contour")
                % "color" is a keyword, but we set it as an explicit argument.
                excludes = [excludes(:); "levels"; "lineSpec"; "color"];
           %elseif strcmp(comp, "heatmap")
           %    excludes = [excludes(:)];
            elseif strcmp(comp, "histfit")
                excludes = [excludes(:); "dist"; "nbins"];
            elseif strcmp(comp, "histogram")
                excludes = [excludes(:); "edges"; "nbins"];
            elseif strcmp(comp, "histogram2")
                excludes = [excludes(:); "nbins"; "xedges"; "yedges"];
            elseif strcmp(comp, "legend")
                excludes = [excludes(:); "labels"];
            elseif strcmp(comp, "plot") || strcmp(comp, "plot3")
                excludes = [excludes(:); "size"]; % color;
            elseif strcmp(comp, "scatter") || strcmp(comp, "scatter3")
                excludes = [excludes(:); "marker"; "filled"; "color"; "size"];
            elseif strcmp(comp, "surface")
                excludes = [excludes(:); "color"; "size"];
            elseif strcmp(comp, "title")
                excludes = [excludes(:); "titletext"; "subtitletext"];
            elseif strcmp(comp, "xlabel") || strcmp(comp, "ylabel") || strcmp(comp, "zlabel")
                excludes = [excludes(:); "txt"];
            elseif ~strcmp(comp, "heatmap")
                disp("comp");
                disp( comp );
                error   ( newline ...
                        + "Internal library error: Unrecognized MATLAB function name" + newline ...
                        + "as ``comp`` argument of object of class ``pm.vis.Axes``." + newline ...
                        + newline ...
                        );
            end
            unique = true;
            onlyfull = true;
            hash = pm.matlab.hashmap.struct2hash(self.(comp), excludes, unique, onlyfull);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end