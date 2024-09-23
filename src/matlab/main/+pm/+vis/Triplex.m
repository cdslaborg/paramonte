%>  \brief
%>  This is the base class for generating instances of
%>  figures containing a square symmetric tiling of subplots.<br>
%>
%>  \details
%>  This class generates figures containing three types of plots in the
%>  upper-triangle, lower-triangle, and diagonal subplots of the figure.<br>
%>  For more information, see the documentation of the class constructor [pm.vis.Triplex.Triplex](@ref Triplex::Triplex).<br>
%>
%>  \note
%>  This class is not meant to be used directly by the end users.<br>
%>  Instead, use the subclasses of this parent class.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.vis.figure.Figure](@ref Figure).<br>
%>
%>  \todo
%>  \phigh
%>  The ``hideShow()`` method of this class still has unresolved bugs that lead to missing X and Y axes labels and tick marks.<br>
%>  Consequently, the ``show()`` and ``hide()`` class methods are also unavailable until the issues are resolved.<br>
%>  For now, all of these functionalities are commented out.<br>
%>  A simple bypass is to generate tilings with prior contemplation about
%>  its design so that post-visualization modifications are unnecessary.<br>
%>
%>  \final
%>
%>  \author
%>  \AmirShahmoradi, August 31 2024, 6:40 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
classdef Triplex < pm.vis.figure.Figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``layout``
        %>
        %>  The scalar object of class [pm.vis.TilingLayout](@ref TilingLayout)
        %>  containing the default layout of the Triplex plot.<br>
        %>  The information includes the subplots and colorbar positions.<br>
        %>  This information is applied only to figure
        %>  components whose positions are unset.<br>
        %>
        layout = [];
        %>
        %>  ``diago``
        %>
        %>  The scalar object of superclass [pm.vis.Subplot](@ref Subplot)
        %>  representing the template of the diagonal subplots to display.<br>
        %>  Note that only the visualization properties of the template are used.<br>
        %>  The data properties of the template are set by the ``make()`` method.<br>
        %>  (**optional**. The default value is [pm.vis.SubplotHistogram](@ref SubplotHistogram).)
        %>
        diago = [];
        %>
        %>  ``lower``
        %>
        %>  The scalar object of superclass [pm.vis.Subplot](@ref Subplot)
        %>  representing the template of the lower-triangle subplots to display.<br>
        %>  The data properties of the template are set by the ``make()`` method.<br>
        %>  (**optional**. The default value is [pm.vis.SubplotContour](@ref SubplotContour).)
        %>
        lower = [];
        %>
        %>  ``upper``
        %>
        %>  The scalar object of superclass [pm.vis.Subplot](@ref Subplot)
        %>  representing the template of the upper-triangle subplots to display.<br>
        %>  The data properties of the template are set by the ``make()`` method.<br>
        %>  (**optional**. The default value is [pm.vis.SubplotLineScatter](@ref SubplotLineScatter).)
        %>
        upper = [];
        %>
        %>  ``tile``
        %>
        %>  The MATLAB cell matrix containing objects of superclass [pm.vis.Subplot](@ref Subplot)
        %>  each of which represents one tile (subplot) axes to display in the figure.<br>
        %>
        tile = cell(0, 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public, Hidden)
        %>
        %>  ``diagYTick``
        %>
        %>  The MATLAB ``cell`` vector of length `nrow``, containing the diagonal y-tick labels.<br>
        %>
        diagYTick = struct();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Construct and return an object of class [pm.vis.Triplex](@ref Triplex).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.Triplex](@ref Triplex).<br>
        %>  The input dataset to each plot section is assumed to be common among all plots (though not necessarily).<br>
        %>  As such, all axes labels and tick marks of all subplots are dropped by default,
        %>  except for subplots in the the left and bottom boundaries of the figure.<br>
        %>
        %>  If any of the three input subplot types has colorbar, it is disabled for individual subplots.<br>
        %>  Instead, universal colorbar(s) will be added to the right and top sides of the Triplex plot.<br>
        %>
        %>  \param[in]  lower       :   The input scalar MATLAB object of superclass [pm.vis.Subplot](@ref Subplot) containing information for
        %>                              the kind the subplot that must be added to the lower-triangular elements of the Triplex matrix of subplots.<br>
        %>                              (**optional**, default = ``[]``. If empty, no lower-triangular subplots will be visualized.)
        %>  \param[in]  diago       :   The input scalar MATLAB object of superclass [pm.vis.Subplot](@ref Subplot) containing information for
        %>                              the kind the subplot that must be added to the diagonal elements of the Triplex matrix of subplots.<br>
        %>                              (**optional**, default = ``[]``. If empty, no diagonal subplots will be visualized.)
        %>  \param[in]  upper       :   The input scalar MATLAB object of superclass [pm.vis.Subplot](@ref Subplot) containing information for
        %>                              the kind the subplot that must be added to the upper-triangular elements of the Triplex matrix of subplots.<br>
        %>                              (**optional**, default = ``[]``. If empty, no upper-triangular subplots will be visualized.)
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.Triplex](@ref Triplex).<br>
        %>
        %>  \interface{Triplex}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(diago, [], [], varargin);
        %>      g = pm.vis.Triplex([], [], upper, varargin);
        %>      g = pm.vis.Triplex([], lower, [], varargin);
        %>
        %>      g = pm.vis.Triplex(diago, lower, [], varargin);
        %>      g = pm.vis.Triplex(diago, [], upper, varargin);
        %>      g = pm.vis.Triplex([], lower, upper, varargin);
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  The data columns to be plotted in each Triplex section must be consistent with other sections.<br>
        %>  Failure to ensure this will lead to undefined behavior and possibly a runtime error.<br>
        %>
        %>  \note
        %>  See the list of class attributes below,
        %>  also those of the superclass [pm.vis.figure.Figure](@ref Figure).<br>
        %>
        %>  \example{Triplex}
        %>  \include{lineno} example/vis/Triplex/main.m
        %>  \vis{Triplex}
        %>  <br><br>
        %>  \image html example/vis/Triplex/Triplex.1.png width=700
        %>  <br><br>
        %>  \image html example/vis/Triplex/Triplex.2.png width=700
        %>  <br><br>
        %>  \image html example/vis/Triplex/Triplex.3.png width=700
        %>  <br><br>
        %>  \image html example/vis/Triplex/Triplex.4.png width=700
        %>  <br><br>
        %>  \image html example/vis/Triplex/Triplex.5.png width=700
        %>  <br><br>
        %>  \image html example/vis/Triplex/Triplex.6.png width=700
        %>
        %>  \final{Triplex}
        %>
        %>  \author
        %>  \AmirShahmoradi, August 31 2024, 6:40 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function self = Triplex(lower, diago, upper, varargin)

            varargin = {"lower", lower, "diago", diago, "upper", upper, varargin{:}};
            self = self@pm.vis.figure.Figure(varargin{:});

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the figure to the original default settings.<br>
        %>
        %>  \details
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.<br>
        %>
        %>  \param[in]  self        :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>      g.reset(varargin); % reset all object properties to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:25 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function reset(self, varargin)

            reset@pm.vis.figure.Figure(self, varargin{:});

            %%%%
            %%%% Ensure equal number of data columns are specified.
            %%%%

            nrow = [];
            ncol = [];
            shape = struct("lower", [], "upper", []);
            for comp = ["lower", "upper"]
                if isa(self.(comp), "pm.vis.SubplotEllipsoid3")
                    shape.(comp) = [numel(self.(comp).dimy), numel(self.(comp).dimx)];
                elseif isa(self.(comp), "pm.vis.Subplot") && isprop(self.(comp), "coly")
                    shape.(comp) = [numel(self.(comp).coly), numel(self.(comp).colx)];
                elseif ~isempty(self.(comp))
                    help("pm.vis.Triplex");
                    disp("class(self.(comp))");
                    disp( class(self.(comp)) );
                    error   ( newline ...
                            + "The specified upper/lower triangular component must be a 2D X-Y" + newline ...
                            + "subplot object of superclass ``pm.vis.SubplotEllipsoid3`` or ``pm.vis.Subplot``." + newline ...
                            + "For more information, see the class documentation displayed above." + newline ...
                            + newline ...
                            );
                end
                if ~isempty(self.(comp))
                    nrow = shape.(comp)(1);
                    ncol = shape.(comp)(2);
                end
            end

            if ~isempty(shape.lower) && ~isempty(shape.upper) % lower and upper both present.
                if ~all(shape.lower == shape.upper)
                    help("pm.vis.Triplex");
                    disp("shape.lower");
                    disp( shape.lower );
                    disp("shape.upper");
                    disp( shape.upper );
                    error   ( newline ...
                            + "The corresponding number of x and y columns in the specified input" + newline ...
                            + "``lower``, ``diago``, and ``upper`` subplot templates must be equal." + newline ...
                            + "For more information, see the class documentation displayed above." + newline ...
                            + newline ...
                            );
                end
            elseif ~isempty(nrow) && ~isempty(ncol) && isa(self.diago, "pm.vis.Subplot") % either lower or upper present.
                if  numel(self.diago.colx) ~= min(nrow, ncol)
                    help("pm.vis.Triplex");
                    disp("[nrow, ncol]");
                    disp( [nrow, ncol] );
                    disp("numel(self.diago.colx)");
                    disp( numel(self.diago.colx) );
                    error   ( newline ...
                            + "The condition ``numel(self.diago.colx) == min(nrow, ncol)`` must hold." + newline ...
                            + "For more information, see the class documentation displayed above." + newline ...
                            + newline ...
                            );
                end
            elseif isa(self.diago, "pm.vis.Subplot") % neither lower nor upper present.
                nrow = numel(self.diago.colx);
                ncol = numel(self.diago.colx);
            elseif ~isempty(self.diago)
                help("pm.vis.Triplex");
                error   ( newline ...
                        + "The input ``diago`` must be 1D plotting object of type number of x and y columns in the specified input" + newline ...
                        + "``lower``, ``diago``, and ``upper`` subplot templates must be equal." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            if ~isempty(nrow) && ~isempty(ncol)
                self.layout = pm.vis.TilingLayout(nrow, ncol);
                self.layout.cbarv.enabled = [];
                self.layout.cbarh.enabled = [];
            else
                help("pm.vis.Triplex");
                disp("class(self.lower)");
                disp( class(self.lower) );
                disp("class(self.diago)");
                disp( class(self.diago) );
                disp("class(self.upper)");
                disp( class(self.upper) );
                error   ( newline ...
                        + "The Triplex class constructor could not determine" + newline ...
                        + "the tiling shape from the specified input arguments." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the make() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.premake(varargin{:}); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Preset the tiling settings before making it.<br>
        %>
        %>  \warning
        %>  This method causes side-effects by manipulating
        %>  the existing attributes of the object.<br>
        %>
        %>  \param[in]  self        :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>      g.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:28 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function premake(self, varargin)

            premake@pm.vis.figure.Figure(self, varargin{:});

            %%%%
            %%%% Reposition figure.
            %%%%

            if  isempty(self.figure.position)
                self.figure.position = self.layout.position;
            end

            %%%%
            %%%% Update the fontsize settings.
            %%%%

            for comp = ["lower", "diago", "upper"]
                if  isa(self.(comp), "pm.vis.Subplot")
                    for prop = ["axes", "xlabel", "ylabel", "zlabel", "colorbar"]
                        if  isprop(self.(comp), prop)
                            if ~isfield(self.(comp).(prop), "fontSize") || isempty(self.(comp).(prop).fontSize)
                                self.(comp).(prop).fontSize = 12;
                            end
                        end
                    end
                end
            end

            %%%%
            %%%% Update the tiling layout to match the latest settings.
            %%%%

            if  isempty(self.layout.cbarv.enabled)
                self.layout.cbarv.enabled = ~isempty(self.lower);
                if  self.layout.cbarv.enabled && isempty(self.lower.cenabled)
                    self.lower.premake();
                    self.layout.cbarv.enabled = self.lower.cenabled && ~self.lower.colorbar.enabled;
                end
                if ~self.layout.cbarv.enabled && ~isempty(self.upper)
                    if  isempty(self.upper.cenabled)
                        self.upper.premake();
                    end
                    self.layout.cbarv.enabled = self.upper.cenabled && ~self.upper.colorbar.enabled;
                end
            end
            if  isempty(self.layout.cbarh.enabled)
                %%%% Activate cbarh only after cbarv is activated.
                self.layout.cbarh.enabled = self.layout.cbarv.enabled && ~isempty(self.lower) && ~isempty(self.upper);
                self.layout.cbarh.enabled = self.layout.cbarh.enabled && self.lower.cenabled && ~self.lower.colorbar.enabled;
                if  self.layout.cbarh.enabled
                    if  isempty(self.upper.cenabled)
                        self.upper.premake();
                    end
                    self.layout.cbarh.enabled = self.upper.cenabled && ~self.upper.colorbar.enabled;
                end
            end
            self.layout.update();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the figure settings and specifications,
        %>  make the figure and the subplots, and return nothing.<br>
        %>
        %>  \details
        %>  The subplots are made by calling their ``make()`` methods.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[in]  self        :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>      g.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function make(self, varargin)

            make@pm.vis.figure.Figure(self, varargin{:});

            self.fout.axes = axes   ( "Parent", self.fout.figure ...
                                    , "position", self.layout.tiling.position ...
                                    , "xlim", [0, 1], "ylim", [0, 1] ...
                                    , "visible", "off" ...
                                    , "color", "none" ...
                                    );

            %%%%
            %%%% Define the subplot cell matrix.
            %%%%

            iplt = 0;
            timer = pm.timing.Timer();
            spinner = pm.timing.Spinner();
            self.tile = cell(self.layout.tiling.nrow, self.layout.tiling.ncol);
            self.diagYTick = cell(min(self.layout.tiling.nrow, self.layout.tiling.ncol), 1);
            for icol = 1 : self.layout.tiling.ncol
                for irow = 1 : self.layout.tiling.nrow

                    iplt = iplt + 1;
                    spinner.spin(iplt / numel(self.tile));

                    %%%%
                    %%%% Copy data from template to the current tile.
                    %%%%

                    comp = [];
                    if  icol < irow && ~isempty(self.lower)
                        comp = "lower";
                    elseif icol == irow && ~isempty(self.diago)
                        comp = "diago";
                    elseif icol > irow && ~isempty(self.upper)
                        comp = "upper";
                    end

                    if ~isempty(comp)

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%% The byte stream approach leads to serious problems with
                        %%%% figures when generated from within sampler components.
                        %self.tile{irow, icol} = getArrayFromByteStream(getByteStreamFromArray(self.(comp)));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        self.tile{irow, icol} = pm.matlab.copy(self.(comp), eval(string(class(self.(comp))) + "()"));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%
                        %%%% Set the data columns to plot in the current tile.
                        %%%%

                        if  isa(self.tile{irow, icol}, "pm.vis.SubplotEllipse3")
                            self.tile{irow, icol}.dimx = self.tile{irow, icol}.dimx(icol);
                            self.tile{irow, icol}.dimy = self.tile{irow, icol}.dimy(irow);
                        elseif isa(self.tile{irow, icol}, "pm.vis.Subplot")
                            if  isprop(self.tile{irow, icol}, "colx")
                                self.tile{irow, icol}.colx = self.tile{irow, icol}.colx(icol);
                            end
                            if  isprop(self.tile{irow, icol}, "coly")
                                self.tile{irow, icol}.coly = self.tile{irow, icol}.coly(irow);
                            end
                        elseif ~isempty(self.tile{irow, icol})
                            error   ( newline ...
                                    + "Unrecognized/unsupported class: " + string(class(self.tile{irow, icol})) + newline ...
                                    + "For more information, see the documentation displayed above." + newline ...
                                    + newline ...
                                    );
                        end

                        if ~isempty(self.tile{irow, icol})

                            %%%%
                            %%%% Set the current tile position.
                            %%%%

                            self.tile{irow, icol}.axes.position = self.layout.tiling.tile.position(:, irow, icol);
                            self.tile{irow, icol}.axes.visible = "on";
                            self.tile{irow, icol}.axes.box = "on";

                            %%%%
                            %%%% Create the current tile axes.
                            %%%%

                            axes_kws = self.tile{irow, icol}.comp2hash("axes");
                            axes(axes_kws{:});

                            %%%%
                            %%%% Make the subplot.
                            %%%%

                            self.tile{irow, icol}.make();

                            %%%%
                            %%%% Set the diagonal YTicks
                            %%%%

                            if  irow == icol
                                xlower = self.tile{irow, icol}.fout.axes.XLim(1);
                                xupper = self.tile{irow, icol}.fout.axes.XLim(2);
                                xrange = xupper - xlower;
                                ylower = self.tile{irow, icol}.fout.axes.YLim(1);
                                yupper = self.tile{irow, icol}.fout.axes.YLim(2);
                                yrange = yupper - ylower;
                                self.diagYTick{irow}.original = struct();
                                self.diagYTick{irow}.modified = struct();
                                self.diagYTick{irow}.original.values = self.tile{irow, icol}.fout.axes.YTick;
                                self.diagYTick{irow}.original.labels = self.tile{irow, icol}.fout.axes.YTickLabel;
                                self.diagYTick{irow}.modified.values = (self.tile{irow, icol}.fout.axes.XTick - xlower) * yrange / xrange + ylower;
                                self.diagYTick{irow}.modified.labels = self.tile{irow, icol}.fout.axes.XTickLabel;
                            end

                        end

                    end

                end
            end

            %%%%
            %%%% Reset the marginal axes labels.
            %%%%

            self.hideShowAxesLabels();

            %%%%
            %%%% Add colorbars.
            %%%%

            shown = struct();
            shown.cbarh = false;
            shown.cbarv = false;
            indexlist = [ self.layout.tiling.nrow, 1 ... lower
                        ; 1, self.layout.tiling.ncol ... upper
                        ; self.layout.tiling.nrow, self.layout.tiling.ncol ... diago
                        ]';
            location = struct();
            location.cbarv = "east";
            location.cbarh = "north";
            for icomp = 1 : size(indexlist, 2)
                irow = indexlist(1, icomp);
                icol = indexlist(2, icomp);
                if ~isempty(self.tile{irow, icol}) && isprop(self.tile{irow, icol}, "cenabled") && self.tile{irow, icol}.cenabled
                    for cbar = ["cbarv", "cbarh"]
                        if ~shown.(cbar) && self.layout.(cbar).enabled

                            %%%%
                            %%%% Make a dummy plot with colorbar from the corresponding tile in the main axes.
                            %%%%

                            kws = {"position", self.layout.tiling.position, "xlim", [0, 1], "ylim", [0, 1], "visible", "off", "color", "none"};
                            dumaxes = axes("Parent", self.fout.figure, kws{:});

                            set(self.fout.figure, 'CurrentAxes', dumaxes);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%% The byte stream approach leads to serious problems with
                            %%%% figures when generated from within sampler components.
                            %dumplot = getArrayFromByteStream(getByteStreamFromArray(self.tile{irow, icol}));
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            dumplot = pm.matlab.copy(self.tile{irow, icol}, eval(string(class(self.tile{irow, icol})) + "()"), [], "fout");
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            dumplot.colorbar.position = self.layout.(cbar).position;
                            dumplot.colorbar.location = location.(cbar);
                            dumplot.colorbar.enabled = true;
                            if  isprop(dumplot, "axes") % Heatmap plots do not have axes property.
                                dumplot.make("axes", kws)
                            else
                                dumplot.make("axes")
                            end
                            dumaxes.Visible = "off";
                            set([dumaxes; dumaxes.Children], 'Visible', 'off');

                            %kws = self.tile{irow, icol}.comp2hash("colorbar");
                            %self.fout.(cbar) = colorbar(dumaxes, kws{:}, "position", self.layout.(cbar).position, "location", location.(cbar));
                            %dfcopy = self.tile{irow, icol}.df.copy();
                            %if ~isprop(self.tile{irow, icol}, "colc")
                            %    colnamc = ["Density"];
                            %elseif ~isempty(self.tile{irow, icol}.colc)
                            %    [~, colnamc] = pm.str.locname(dfcopy.Properties.VariableNames, self.tile{irow, icol}.colc);
                            %else
                            %    colnamc = "Row Index";
                            %end
                            %xlabel(self.fout.(cbar), colnamc(1));
                            %cmap = colormap(self.tile{irow, icol}.fout.axes);
                            %if ~isfield(self.fout, "colormap")
                            %   self.fout.colormap = struct();
                            %end
                            %self.fout.colormap.(cbar) = colormap(self.fout.(cbar), cmap);
                            shown.(cbar) = true;
                            break;

                        end
                    end
                    if  shown.cbarh && shown.cbarv
                        break;
                    end
                end
            end

            disp("done in " + sprintf("%.6f", string(timer.toc())) + " seconds.");

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Rotate the axes labels of the subplots of the Triplex, and return nothing.<br>
        %>
        %>  \details
        %>  All axes labels in the tiling will be impacted by this routine,
        %>  even though only a selected subset might be visible.<br>
        %>
        %>  \param[in]  self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>  \param[in]  degx    :   The input scalar MATLAB positive whole-number, representing the amount of
        %>                          rotation to be applied to the x-axis labels with respect to the horizontal line.<br>
        %>                          If it is set to empty ``[]``, the axis label orientation will remain intact.<br>
        %>                          (**optional**, default = ``45``. It must be present if and only if ``degy`` is also present.)
        %>  \param[in]  degy    :   The input scalar MATLAB positive whole-number, representing the amount of
        %>                          rotation to be applied to the y-axis labels with respect to the horizontal line.<br>
        %>                          If it is set to empty ``[]``, the axis label orientation will remain intact.<br>
        %>                          (**optional**, default = ``45``. It must be present if and only if ``degx`` is also present.)
        %>
        %>  \interface{rotateAxesLabels}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>      g.make(varargin);
        %>
        %>      g.rotateAxesLabels(); % rotate all axes labels by 45 degrees.
        %>      g.rotateAxesLabels(degx, []); % rotate x-axis labels by ``degx`` degrees.
        %>      g.rotateAxesLabels([], degy); % rotate y-axis labels by ``degy`` degrees.
        %>      g.rotateAxesLabels(degx, degy); % rotate x-axis and y-axis labels by ``degx`` and ``degy`` degrees.
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  This method has side-effects by manipulating the existing attributes of the parent object.<br>
        %>
        %>  \final{rotateAxesLabels}
        %>
        %>  \author
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function rotateAxesLabels(self, degx, degy)
            if  nargin < 2
                degx = 45;
                degy = 45;
            elseif nargin ~= 3
                help("pm.vis.Triplex.rotateAxesLabels");
                error   ( newline ...
                        + "The ``rotateAxesLabels()`` method takes either no or exactly two input arguments." + newline ...
                        + "If no argument is passed, axes labels will be rotated 45 degrees." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            for icol = 1 : size(self.tile, 2)
                for irow = 1 : size(self.tile, 2)
                    self.tile{irow,icol}.fout.axes.XLabel.Rotation = degx;
                    self.tile{irow,icol}.fout.axes.YLabel.Rotation = degy;
                    self.tile{irow,icol}.fout.axes.XLabel.VerticalAlignment = "top";
                    self.tile{irow,icol}.fout.axes.YLabel.VerticalAlignment = "middle";
                    self.tile{irow,icol}.fout.axes.XLabel.HorizontalAlignment = "right";
                    self.tile{irow,icol}.fout.axes.YLabel.HorizontalAlignment = "right";
                end
            end
        end % rotateAxesLabels

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%>  \brief
        %%%>  Hide the requested section(s) of the Triplex, and return nothing.<br>
        %%%>
        %%%>  \details
        %%%>  Note that the main axes colorbar associated with the specified component is also automatically hidden.<br>
        %%%>
        %%%>  \param[in]  self        :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %%%>  \param[in]  complist    :   The input vector of MATLAB strings each of which can be any of the following:
        %%%>                              <ol>
        %%%>                                  <li>    ``"lower"``     :   corresponding to the lower triangle of the Triplex plot.
        %%%>                                  <li>    ``"diago"``     :   corresponding to the diagonal tiles in the Triplex plot.
        %%%>                                  <li>    ``"upper"``     :   corresponding to the upper triangle of the Triplex plot.
        %%%>                              </ol>
        %%%>                              (**optional**, default = ``["lower", "diago", "upper"]``.)<br>
        %%%>
        %%%>  \interface{hide}
        %%%>  \code{.m}
        %%%>
        %%%>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %%%>      g.make(varargin);
        %%%>      g.hide(complist);
        %%%>      g.hide([]);
        %%%>
        %%%>  \endcode
        %%%>
        %%%>  \warning
        %%%>  This method has side-effects by manipulating the existing attributes of the parent object.<br>
        %%%>
        %%%>  \example{hide}
        %%%>  \code{.m}
        %%%>
        %%%>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %%%>      g.make(varargin);
        %%%>      g.hide() % hide all parts of the Triplex plot and the associated colorbars in the Triplex plot.
        %%%>      g.hide(["diago", "upper"]) % hide the diagonal subplots and the associated colorbars in the Triplex plot.
        %%%>      g.hide("lower") % hide the lower triangle subplots of the Triplex plot and its associated colorbar in the Triplex plot.
        %%%>
        %%%>  \endcode
        %%%>
        %%%>  \final{hide}
        %%%>
        %%%>  \author
        %%%>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %%function hide(self, complist)
        %%    if  nargin < 2
        %%        complist = [];
        %%    end
        %%    self.hideShow("hide", complist);
        %%end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        %%%>  \brief
        %%%>  Show the requested section(s) of the Triplex, and return nothing.<br>
        %%%>
        %%%>  \details
        %%%>  Note that the main axes colorbar associated with the specified component is also automatically shown.<br>
        %%%>
        %%%>  \param[in]  complist    :   The input vector of MATLAB strings each of which can be any of the following:
        %%%>                              <ol>
        %%%>                                  <li>    ``"lower"``     :   corresponding to the lower triangle of the Triplex plot.
        %%%>                                  <li>    ``"diago"``     :   corresponding to the diagonal tiles in the Triplex plot.
        %%%>                                  <li>    ``"upper"``     :   corresponding to the upper triangle of the Triplex plot.
        %%%>                              </ol>
        %%%>                              (**optional**, default = ``["lower", "diago", "upper"]``.)<br>
        %%%>
        %%%>  \interface{show}
        %%%>  \code{.m}
        %%%>
        %%%>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %%%>      g.make(varargin);
        %%%>      g.show(complist);
        %%%>
        %%%>  \endcode
        %%%>
        %%%>  \warning
        %%%>  This method has side-effects by manipulating the existing attributes of the parent object.<br>
        %%%>
        %%%>  \example{show}
        %%%>  \code{.m}
        %%%>
        %%%>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %%%>      g.make(varargin);
        %%%>      g.show() % show all parts of the Triplex plot and the associated colorbars in the Triplex plot.
        %%%>      g.show(["diago", "upper"]) % show the diagonal subplots and the associated colorbars in the Triplex plot.
        %%%>      g.show("lower") % show the lower triangle subplots of the Triplex plot and its associated colorbar in the Triplex plot.
        %%%>
        %%%>  \endcode
        %%%>
        %%%>  \final{show}
        %%%>
        %%%>  \author
        %%%>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %%function show(self, complist)
        %%    if  nargin < 2
        %%        complist = [];
        %%    end
        %%    self.hideShow("show", complist);
        %%end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%>  \brief
        %%%>  Hide/Show the requested section(s) of the Triplex, and return nothing.<br>
        %%%>
        %%%>  \details
        %%%>  This is an internal (hidden) method of the class and inaccessible to end users.<br>
        %%%>  Note that the main axes colorbar associated with the specified component is also automatically impacted.<br>
        %%%>
        %%%>  \param[in]  self        :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %%%>  \param[in]  action      :   The input scalar MATLAB string that can be any of the following:<br>
        %%%>                              <ol>
        %%%>                                  <li>    ``"hide"``      :   hide the specified components of the Triplex plot in the input argument ``complist``.
        %%%>                                  <li>    ``"show"``      :   show the specified components of the Triplex plot in the input argument ``complist``.
        %%%>                              </ol>
        %%%>                              (**optional**, default = ``["lower", "diago", "upper"]``.)<br>
        %%%>  \param[in]  complist    :   The input vector of MATLAB strings each of which can be any of the following:
        %%%>                              <ol>
        %%%>                                  <li>    ``"lower"``     :   corresponding to the lower triangle of the Triplex plot.
        %%%>                                  <li>    ``"diago"``     :   corresponding to the diagonal tiles in the Triplex plot.
        %%%>                                  <li>    ``"upper"``     :   corresponding to the upper triangle of the Triplex plot.
        %%%>                              </ol>
        %%%>                              (**optional**, default = ``["lower", "diago", "upper"]``.)<br>
        %%%>
        %%%>  \interface{show}
        %%%>  \code{.m}
        %%%>
        %%%>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %%%>      g.make(varargin);
        %%%>      g.hideShow(action, []);
        %%%>      g.hideShow(action, complist);
        %%%>
        %%%>  \endcode
        %%%>
        %%%>  \warning
        %%%>  This method has side-effects by manipulating the existing attributes of the parent object.<br>
        %%%>
        %%%>  \final{show}
        %%%>
        %%%>  \author
        %%%>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %%function hideShow(self, action, complist)
        %%
        %%    action = string(action);
        %%    showit = strcmpi(action, "show");
        %%    hideit = strcmpi(action, "hide");
        %%    if  showit
        %%        onoff = "on";
        %%    elseif hideit
        %%        onoff = "off";
        %%    end
        %%
        %%    %%%%
        %%    %%%% Determine which figure parts to act on.
        %%    %%%%
        %%
        %%    if  nargin == 2 || isempty(complist)
        %%        complist = ["lower", "diago", "upper"];
        %%    end
        %%
        %%    %%%%
        %%    %%%% Act on the specific requested parts of the figure.
        %%    %%%%
        %%
        %%    enabled = struct();
        %%    enabled.diago = false;
        %%    enabled.lower = false;
        %%    enabled.upper = false;
        %%    complist = lower(string(complist));
        %%    for comp = complist
        %%        if strcmpi(comp, "diag") || strcmpi(comp, "diago")
        %%            enabled.diago  = true;
        %%        elseif strcmpi(comp, "upper")
        %%            enabled.upper = true;
        %%        elseif strcmpi(comp, "lower")
        %%            enabled.lower = true;
        %%        else
        %%            help("pm.vis.Triplex." + action);
        %%            error   ( newline ...
        %%                    + "Unrecognized input argument pass to the " + action + "() method: " + newline ...
        %%                    + newline ...
        %%                    + comp + newline ...
        %%                    + newline ...
        %%                    + "The " + action + "() method accepts only a list of comma-separated" + newline ...
        %%                    + "string values representing the ""part"" of the figure to " + action + "." + newline ...
        %%                    + "For more information, see the documentation of the method displayed above." + newline ...
        %%                    + newline ...
        %%                    );
        %%        end
        %%    end
        %%
        %%    %%%%
        %%    %%%% Perform action.
        %%    %%%%
        %%
        %%    for icol = 1 : size(self.tile, 2)
        %%        for irow = 1 : size(self.tile, 1)
        %%            if (icol < irow && enabled.lower) || (irow < icol && enabled.upper) || (icol == irow && enabled.diago)
        %%                if ~isempty(self.tile{irow, icol}.axes)
        %%                    self.tile{irow, icol}.fout.axes.Visible = onoff;
        %%                    set(get(self.tile{irow, icol}.fout.axes, "children"), "visible", onoff);
        %%                    %if  icol < irow && self.layout.cbarv.enabled
        %%                    %    if  strcmpi(onoff, "hide")
        %%                    %    end
        %%                    %end
        %%                end
        %%            end
        %%        end
        %%    end
        %%
        %%    if enabled.lower || enabled.upper || enabled.diago
        %%        self.hideShowAxesLabels();
        %%    end
        %%
        %%end % function hideShow

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Hide or show axis labels and ticks depending on the presence of the neighbor subplots.<br>
        %>
        %>  \details
        %>  This is an internal (hidden) method of the class and inaccessible to end users.<br>
        %>
        %>  \param[in]  self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>
        %>  \interface{hideShowAxesLabels}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>      g.hideShowAxesLabels();
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  This method causes side-effects by manipulating the existing attributes of the object.<br>
        %>
        %>  \final{hideShowAxesLabels}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:28 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function hideShowAxesLabels(self)

            %%%%
            %%%% (Un)set the tile axes labels and ticks.
            %%%%

            for icol = 1 : size(self.tile, 2)
                for irow = 1 : size(self.tile, 1)

                    if ~isempty(self.tile{irow, icol})% && strcmpi(self.tile{irow, icol}.fout.axes.XLabel.Visible, "on")

                        zisabled = false;
                        if  irow < size(self.tile, 1)
                            if ~isempty(self.tile{irow + 1, icol}) && strcmpi(self.tile{irow + 1, icol}.fout.axes.XLabel.Visible, "on")
                                self.tile{irow, icol}.fout.axes.XLabel.Visible = "off";
                               %self.tile{irow, icol}.fout.axes.XLabel.String = "";
                                self.tile{irow, icol}.fout.axes.XTickLabel = [];
                                zisabled = true;
                            elseif isempty(self.tile{irow + 1, icol})
                                self.tile{irow, icol}.fout.axes.XTickLabelMode = "auto";
                            %    self.tile{irow, icol}.fout.axes.XLabel.Visible = "on";
                            end
                        end

                        if  icol > 1
                            if ~isempty(self.tile{irow, icol - 1}) && strcmpi(self.tile{irow, icol - 1}.fout.axes.Visible, "on")
                                self.tile{irow, icol}.fout.axes.YLabel.Visible = "off";
                               %self.tile{irow, icol}.fout.axes.YLabel.String = "";
                                self.tile{irow, icol}.fout.axes.YTickLabel = [];
                                zisabled = true;
                            elseif isempty(self.tile{irow, icol - 1})
                                self.tile{irow, icol}.fout.axes.YTickLabelMode = "auto";
                            %    self.tile{irow, icol}.fout.axes.YLabel.Visible = "on";
                            end
                        end

                        if  zisabled
                            self.tile{irow, icol}.fout.axes.ZLabel.Visible = "off";
                           %self.tile{irow, icol}.fout.axes.ZLabel.String = "";
                            self.tile{irow, icol}.fout.axes.ZTickLabel = [];
                        else
                            self.tile{irow, icol}.fout.axes.ZTickLabelMode = "auto";
                        %    self.tile{irow, icol}.fout.axes.ZLabel.Visible = "off";
                        end

                        if  icol ~= 1
                            self.tile{irow, icol}.fout.axes.YLabel.Rotation = 45;
                            self.tile{irow, icol}.fout.axes.YLabel.VerticalAlignment = "middle";
                            self.tile{irow, icol}.fout.axes.YLabel.HorizontalAlignment = "right";
                            self.tile{irow, icol}.fout.axes.YLabel.FontSize = self.tile{irow, icol}.fout.axes.FontSize;
                        end

                        if  irow ~= size(self.tile, 1)
                            self.tile{irow, icol}.fout.axes.XLabel.Rotation = 45;
                            self.tile{irow, icol}.fout.axes.XLabel.VerticalAlignment = "top";
                            self.tile{irow, icol}.fout.axes.XLabel.HorizontalAlignment = "right";
                            self.tile{irow, icol}.fout.axes.XLabel.FontSize = self.tile{irow, icol}.fout.axes.FontSize;
                        end

                    end

                end
            end

            self.setDiagoTicksY();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Hide or show axis labels and ticks depending on the presence of the neighbor subplots.<br>
        %>
        %>  \details
        %>  This is an internal (hidden) method of the class and inaccessible to end users.<br>
        %>
        %>  \param[in]  self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>
        %>  \interface{hideShowAxesLabels}
        %>  \code{.m}
        %>
        %>      g = pm.vis.Triplex(lower, diago, upper, varargin);
        %>      g.hideShowAxesLabels();
        %>
        %>  \endcode
        %>
        %>  \warning
        %>  This method causes side-effects by manipulating the existing attributes of the object.<br>
        %>
        %>  \final{hideShowAxesLabels}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:28 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function setDiagoTicksY(self)

            %%%%
            %%%% reset the ytick labels of the diagonal subplots to x-axis ticks when there are subplots to their right.
            %%%%

            for iell = 1 : min(size(self.tile))

                if ~isempty(self.tile{iell, iell}) && strcmp(self.tile{iell, iell}.fout.axes.Visible, "on")

                    %%%%
                    %%%% Do NOT change the order of conditions in the if blocks. short-circuit is at work.
                    %%%%

                    if (iell == size(self.tile, 2) && (isempty(self.tile{iell, iell - 1}) || strcmp(self.tile{iell, iell - 1}.fout.axes.Visible, "off"))) ...
                    || (iell == 1 && (isempty(self.tile{iell, iell + 1}) || strcmp(self.tile{iell, iell + 1}.fout.axes.Visible, "off"))) ...
                    || (iell ~= 1 && iell ~= size(self.tile, 2) ...
                        && ...
                        ( isempty(self.tile{iell, iell + 1}) || strcmp(self.tile{iell, iell + 1}.fout.axes.Visible, "off")) ...
                        && ...
                        ( isempty(self.tile{iell, iell - 1}) || strcmp(self.tile{iell, iell - 1}.fout.axes.Visible, "off")) ...
                        )

                        %%%%
                        %%%% Let the histograms have their own y-ranges.
                        %%%%

                        %set(self.tile{iell, iell}.fout.axes, "YTickMode", "auto", "YTickLabelMode", "auto");

                    elseif iell < size(self.tile, 2)

                        %%%%
                        %%%% Set the y-range and ticklabels according to the needs of the neighbor plots.
                        %%%%

                        %%%% ``self.tile{iell, iell}.fout`` is the current tile.
                        %%%% ``self.tile{inei, iell}.fout`` is the neighbor tile.

                        if  iell == 1 || isempty(self.tile{iell, iell - 1}) || strcmp(self.tile{iell, iell - 1}.fout.axes.Visible, "off")

                            if ~isempty(self.tile{end, iell})
                                inei = size(self.tile, 1);
                            else
                                inei = iell;
                            end
                            yticks(self.tile{iell, iell}.fout.axes, self.diagYTick{iell}.modified.values);
                            yticklabels(self.tile{iell, iell}.fout.axes, self.diagYTick{iell}.modified.labels);

                            %%%%
                            %%%% Set Y-axis label.
                            %%%%

                            failed = ~isfield(self.tile{inei, iell}.fout, "xlabel") || isempty(self.tile{inei, iell}.fout.xlabel.String);
                            if ~failed
                                try
                                    self.tile{iell, iell}.fout.axes.YLabel.String = self.tile{inei, iell}.fout.xlabel.String;
                                catch
                                    failed = true;
                                end
                            end
                            if  failed
                                failed = isempty(self.tile{iell, iell}.fout.xlabel.String);
                                if ~failed
                                    try
                                        self.tile{iell, iell}.fout.axes.YLabel.String = self.tile{inei, iell}.fout.xlabel.String;
                                    catch
                                        failed = true;
                                    end
                                end
                                if  failed
                                    self.tile{iell, iell}.fout.axes.YLabel.String = string(self.tile{iell, iell}.df.Properties.VariableNames(iell));
                                end
                            end

                        end

                    end

                end

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end