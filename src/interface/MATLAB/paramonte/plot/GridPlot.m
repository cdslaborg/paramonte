classdef GridPlot < BasePlot

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = protected, Hidden)
        %dfref = [];
        %isdryrun = [];
        plotTypeList =  [ "line" ...
                        , "hist" ...
                        , "hist2" ...
                        , "histfit" ...
                        , "scatter" ...
                        , "scatter3" ...
                        , "lineScatter" ...
                        , "lineScatter3" ...
                        ];
        plotTypeListLen
        lowerEnabled
        upperEnabled
        diagEnabled
    end

    properties (Access = public)

        columns
        ccolumn
        %line_kws
        %hist_kws
        %hist2_kws
        %histfit_kws
        %scatter_kws
        %scatter3_kws
        %lineScatter_kws
        %lineScatter3_kws
        colorbar_kws
        colormap
        layout

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function reset(self)

            reset@BasePlot(self);
            self.plotTypeListLen = length(self.plotTypeList);

            %self.columns = {};
            self.ccolumn = {};
            for i = 1:length(self.plotTypeList)
                self.layout.subplot.(self.plotTypeList{i}+"_kws") = struct();
            end
            self.colormap = {};
            self.colorbar_kws = {};
            self.layout.plotType.upper = "lineScatter";
            self.layout.plotType.lower = "hist2";
            self.layout.plotType.diag = "hist";
            self.layout.update = @self.updateLayout;

            self.layout.axis.main.margin.bottom = 0.07; % 0.14; 
            self.layout.axis.main.margin.right = 0.10;
            self.layout.axis.main.margin.left = 0.05; %0.1;
            self.layout.axis.main.margin.top = 0.0;
            self.layout.axis.main.origin.x = 0.02;
            self.layout.axis.main.origin.y = 0.00;
            %self.layout.axis.main.height = 1 - self.layout.axis.main.margin.top;
            %self.layout.axis.main.width = 1 - self.layout.axis.main.margin.right;
            self.layout.axis.main.title.content = [];
            self.layout.axis.main.title.fontsize = 13;
            self.layout.axis.main.nrow = length(self.columns);
            self.layout.axis.main.ncol = self.layout.axis.main.nrow;

            self.layout.axis.subplot.interspace = 0.015; % space between subplots
            %self.layout.axis.subplot.height = (1-self.layout.axis.main.margin.top-self.layout.axis.main.margin.bottom)/self.layout.axis.main.ncol - self.layout.axis.subplot.interspace;
            %self.layout.axis.subplot.width  = (1-self.layout.axis.main.margin.left-self.layout.axis.main.margin.right)/self.layout.axis.main.nrow - self.layout.axis.subplot.interspace;
            self.layout.axis.subplot.box = "on";
            self.layout.axis.subplot.grid = "on";
            self.layout.axis.subplot.fontsize = 11;

            %self.layout.colorbar.origin.x = 1 - self.layout.axis.main.margin.right;
            %self.layout.colorbar.origin.y = self.layout.axis.main.margin.bottom;
            %self.layout.colorbar.height = self.layout.axis.main.nrow * ( self.layout.axis.subplot.height + self.layout.axis.subplot.interspace ) - self.layout.axis.subplot.interspace;
            self.layout.colorbar.width = 0.03;
            self.layout.colorbar.fontsize = 12;

            self.layout.figHeight   = 900;
            %self.layout.figWidth    = self.layout.figHeight *   ( 1 ...
            %                                                    - self.layout.axis.main.margin.bottom ...
            %                                                    + self.layout.axis.main.margin.right ...
            %                                                    + self.layout.axis.main.margin.left ...
            %                                                    - self.layout.axis.main.margin.top ...
            %                                                    + self.layout.colorbar.width ...
            %                                                    );

            self.updateLayout();

            %self.target    = Target()
            self.isdryrun = true;
            self.plot();
            self.isdryrun = false;

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = GridPlot(varargin) % expected input arguments: dataFrame, plotType

            self = self@BasePlot(varargin{1});
            try
                self.columns = string(varargin{2});
            catch
                error   ( "GridPlot requires two input arguments: " + newline + newline ...
                        + "    the data Table " + newline ...
                        + "    and a cell array of column names to plot " + newline + newline ...
                        + "The column names is required in order to set the layout of the Grid plot." ...
                        );
            end
            self.reset()
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function plot(self,varargin)
        %   Generate a line plot from the selected columns of the object's dataframe.
        %   
        %   Parameters
        %   ----------
        %       None
        %   
        %   Returns
        %   -------
        %       None. However, this method causes side-effects by manipulating 
        %       the existing attributes of the object.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% parse arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseArgs(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% adjust figure position
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.gca_kws = {};
            self.gcf_kws = addKeyVal( "position",   [ 50                    ... origin x
                                                    , 100                   ... origin y
                                                    , self.layout.figWidth  ... width
                                                    , self.layout.figHeight ... height
                                                    ], self.gcf_kws{:} );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.ccolumn)>1
                error( "the ccolumn attribute should be set to a single column name from the input data Table.");
            end

            cEnabled = getVecLen(self.ccolumn) || isa(self.ccolumn,"cell");

            if cEnabled && ~getVecLen(self.colormap)
                self.colormap = "winter";
            end

            % check subplot data columns to plot

            rowsIsPresent = getVecLen(self.rows);
            columnsIsPresent = getVecLen(self.columns);
            ccolumnIsPresent = getVecLen(self.ccolumn);
            for i = 1:self.plotTypeListLen

                plotType_kws = self.plotTypeList(i) + "_kws";

                % check rows presence

                if rowsIsPresent
                    self.layout.subplot.(plotType_kws).rows = self.rows;
                else
                    self.layout.subplot.(plotType_kws).rows = {};
                end

                %self.layout.subplot.(plotType_kws).xcolumns = {};

                % set color properties

                if ~( strcmp(self.plotTypeList(i),"hist") || strcmp(self.plotTypeList(i),"histfit") )

                    %self.layout.subplot.(plotType_kws).ycolumns = {};

                    if ccolumnIsPresent
                        if contains(self.plotTypeList(i),"3")
                            self.layout.subplot.(plotType_kws).zcolumns = self.ccolumn;
                        end
                        if strcmp(self.plotTypeList(i),"hist2")
                            self.layout.subplot.(plotType_kws).histogram2_kws = {"NumBins",[100 100]};
                            self.layout.subplot.(plotType_kws).colormap = {};
                        else
                            self.layout.subplot.(plotType_kws).ccolumns = self.ccolumn;
                            if ~isfield(self.layout.subplot.(plotType_kws),"plot_kws") || isempty(self.layout.subplot.(plotType_kws).plot_kws)
                                self.layout.subplot.(plotType_kws).plot_kws = {};
                            end
                            self.layout.subplot.(plotType_kws).plot_kws = addKeyVal("linewidth", 0.75, self.layout.subplot.(plotType_kws).plot_kws{:});
                            self.layout.subplot.(plotType_kws).plot_kws = addKeyVal("color", uint8([200 200 200]), self.layout.subplot.(plotType_kws).plot_kws{:});
                            self.layout.subplot.(plotType_kws).surface_kws = [];
                            self.layout.subplot.(plotType_kws).colormap = {}; % self.colormap
                        end
                        self.layout.subplot.(plotType_kws).colorbar_kws = [];
                    end

                end

                % check gcf_kws/box/grid presence

                self.layout.subplot.(plotType_kws).legend_kws = [];
                self.layout.subplot.(plotType_kws).gcf_kws = [];
                self.layout.subplot.(plotType_kws).gca_kws =    { "box", self.layout.axis.subplot.box ...
                                                                ..., "grid", self.layout.axis.subplot.grid ...
                                                                , "fontsize", self.layout.axis.subplot.fontsize ...
                                                                };

            end

            % check columns presence

            if ~columnsIsPresent
                try
                    self.columns = self.dfref.Properties.VariableNames;
                catch
                    self.columns = [];
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % generate figure and axes if needed

            figEnabled = isa(self.gcf_kws,"string") || isa(self.gcf_kws,"cell");

            if figEnabled
                self.currentFig.gcf = figure(self.gcf_kws{:});
            else
                self.currentFig.gcf = gcf;
            end
            set(0, "CurrentFigure", self.currentFig.gcf);

            % construct plot handles

            self.lowerEnabled = (getVecLen(self.layout.plotType.lower) > 0) || isa(self.layout.plotType.lower,"cell");
            self.upperEnabled = (getVecLen(self.layout.plotType.upper) > 0) || isa(self.layout.plotType.upper,"cell");
            self.diagEnabled = (getVecLen(self.layout.plotType.diag) > 0) || isa(self.layout.plotType.diag,"cell");

            if self.lowerEnabled
                if isempty(self.layout.plotType.lower)
                    self.layout.plotType.lower = "hist2";
                else
                    self.layout.plotType.lower = string(self.layout.plotType.lower);
                end
            end

            if self.upperEnabled
                if isempty(self.layout.plotType.upper)
                    self.layout.plotType.upper = "lineScatter";
                else
                    self.layout.plotType.upper = string(self.layout.plotType.upper);
                end
            end

            if self.diagEnabled
                if isempty(self.layout.plotType.diag)
                    self.layout.plotType.diag = "hist";
                else
                    self.layout.plotType.diag = string(self.layout.plotType.diag);
                end
            end

            % generate axes and subplots
            %hold on;
            self.currentFig.gca = axes  ( "position" ...
                                        , self.getMainAxisPosition()  ...
                                        ..., "Position", mainPosition ...
                                        , "Xlim", [0 1], "Ylim", [0 1] ...
                                        , "Color", "none" ...
                                        , "Parent", self.currentFig.gcf ...
                                        );
            axis(self.currentFig.gca,"off");
            %set(self.currentFig.gca,"visible","off");
            %hold off;


            axisCount = self.layout.axis.main.nrow * self.layout.axis.main.ncol;
            msgSuffix = " out of " + string(axisCount);

            if self.lowerEnabled 
                if contains(self.layout.plotType.lower,"line") || contains(self.layout.plotType.lower,"scatter")
                    lowerPlotClass = "LineScatterPlot";
                else
                    lowerPlotClass = "HistPlot";
                end
            end
            if self.upperEnabled 
                if contains(self.layout.plotType.upper,"line") || contains(self.layout.plotType.upper,"scatter")
                    upperPlotClass = "LineScatterPlot";
                else
                    upperPlotClass = "HistPlot";
                end
            end
            if self.diagEnabled 
                diagPlotClass = "HistPlot"; %self.layout.plotType.diag;
            end

            iaxis = 0;
            self.currentFig.subplotList = cell(self.layout.axis.main.ncol,self.layout.axis.main.nrow);
            for irow = 1:self.layout.axis.main.nrow %:-1:1
                for icol = 1:self.layout.axis.main.ncol%:-1:1

                    jrow = self.layout.axis.main.nrow-irow+1;
                    jcol = self.layout.axis.main.ncol-icol+1;
                    isLower = irow>icol && self.lowerEnabled;
                    isUpper = irow<icol && self.upperEnabled;
                    isDiag = irow==icol && self.diagEnabled;

                    iaxis = iaxis + 1;

                    if isLower || isUpper || isDiag

                        disp( "generating subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix )

                        currentSubplotAxis = axes   ( "position" ...
                                                    , self.getSubplotAxisPosition(jrow,icol) ...
                                                    );

                        hold on;

                        xcolumns = self.columns(jcol);
                        ycolumns = self.columns(jrow);

                        if isLower
                            currentPlotType = self.layout.plotType.lower;
                            currentPlotTypeClass = lowerPlotClass;
                            self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(self.dfref,currentPlotType)");
                            self.currentFig.subplotList{irow,icol}.xcolumns = xcolumns;
                            self.currentFig.subplotList{irow,icol}.ycolumns = ycolumns;
                        elseif isUpper
                            currentPlotType = self.layout.plotType.upper;
                            currentPlotTypeClass = upperPlotClass;
                            self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(self.dfref,currentPlotType)");
                            self.currentFig.subplotList{irow,icol}.xcolumns = xcolumns;
                            self.currentFig.subplotList{irow,icol}.ycolumns = ycolumns;
                        else
                            currentPlotType = self.layout.plotType.diag;
                            currentPlotTypeClass = diagPlotClass;
                            self.currentFig.subplotList{irow,icol} = eval(currentPlotTypeClass+"(self.dfref,currentPlotType)");
                            self.currentFig.subplotList{irow,icol}.xcolumns = xcolumns;
                        end
                        %self.currentFig.subplotList(iaxis).parseArgs(inputArgs{:});

                        plotType_kws = currentPlotType + "_kws";
                        fieldList = fieldnames(self.layout.subplot.(plotType_kws));
                        valueList = struct2cell(self.layout.subplot.(plotType_kws));
                        for i = 1:length(fieldList)
                            self.currentFig.subplotList{irow,icol}.(fieldList{i}) = valueList{i};
                            %inputArgs = { inputArgs{:} , fieldList{i}, valueList(i) };
                        end

                        if strcmp(currentPlotType,"hist2")
                            %colormap();
                            %self.layout.subplot.(plotType_kws).colormap = flipud(cold());
                            %self.layout.subplot.(plotType_kws).colormap = "autumn";
                            %colormap(self.layout.subplot.(plotType_kws).colormap);
                        end

                        self.currentFig.subplotList{irow,icol}.plot()
                        self.currentFig.subplotList{irow,icol}.currentFig.gca = currentSubplotAxis;

                        if irow==self.layout.axis.main.nrow
                            %axisLabelX = get(gca,'XLabel');
                            %set(axisLabelX,'rotation',45,'VerticalAlignment','top','HorizontalAlignment','right');
                        else
                            set(gca,'XLabel',[]);
                        end
                        if icol==1
                            %axisLabelY = get(gca,'YLabel');
                            %set(axisLabelY,'rotation',45,'VerticalAlignment','middle','HorizontalAlignment','right');
                            if irow==1
                                xlower = self.currentFig.subplotList{1,1}.currentFig.gca.XLim(1);
                                xupper = self.currentFig.subplotList{1,1}.currentFig.gca.XLim(2);
                                xrange = xupper - xlower;
                                ylower = self.currentFig.subplotList{1,1}.currentFig.gca.YLim(1);
                                yupper = self.currentFig.subplotList{1,1}.currentFig.gca.YLim(2);
                                yrange = yupper - ylower;
                                newYTicks = (self.currentFig.subplotList{1,1}.currentFig.gca.XTick - xlower) * yrange / xrange + ylower;
                                yticks(newYTicks);
                                yticklabels(self.currentFig.subplotList{1,1}.currentFig.gca.XTickLabel);
                                ylabel(string(self.currentFig.subplotList{1,1}.xcolumns));
                                %self.currentFig.subplotList{1,1}.currentFig.gca.YTick = self.currentFig.subplotList{1,1}.currentFig.gca.XTick;
                                %self.currentFig.subplotList{1,1}.currentFig.gca.YTickLabel = self.currentFig.subplotList{1,1}.currentFig.gca.XTickLabel;
                            end
                        else
                            set(gca,'YLabel',[]);
                        end
                        set(gca,'ZLabel',[]);

                        set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"FontSize",self.layout.axis.subplot.fontsize);
                        if icol ~= 1
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YTickLabels",[]);
                        end
                        if irow ~= self.layout.axis.main.nrow
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XTickLabels",[]);
                        end

                        %set(self.currentFig.axisList(iaxis)(iaxis),'box',self.layout.axis.subplot.box);
                        hold off;

                    else

                        disp( "skipping subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix )
                        self.currentFig.subplotList{irow,icol} = [];

                    end
                    
                    if strcmp(self.layout.axis.subplot.grid,"on")
                        grid on;
                    elseif strcmp(self.layout.axis.subplot.grid,"off")
                        grid off;
                    end

                end
            end

            % set up colorbar

            if cEnabled

                colormap(self.colormap);

                ccolumnValues = self.dfref.(self.ccolumn);
                colorbarLimit = [ min(ccolumnValues) max(ccolumnValues) ]; % Colorbar range
                colorbarLabel = self.ccolumn;

                caxis(colorbarLimit);
                set(self.currentFig.gca,"CLim",colorbarLimit);
                self.currentFig.colorbar = colorbar();
                ylabel(self.currentFig.colorbar,colorbarLabel);
                set ( self.currentFig.colorbar ...
                    , "position" ...
                    , self.getColorbarAxisPosition() ...
                    , "FontSize", self.layout.colorbar.fontsize ...
                    ) ;

            end

        end % function plot

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function hide(self,varargin)
            self.hideShow(varargin{:},"hide")
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function show(self,varargin)
            self.hideShow(varargin{:},"show")
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function rotateAxisLabels(self,varargin)

            if nargin==1
                degreeX = 45;
                degreeY = 45;
            elseif nargin==2
                degreeX = varargin{1};
                degreeY = varargin{1};
            elseif nargin==3
                degreeX = varargin{1};
                degreeY = varargin{2};
            else
                error   ( newline ...
                        + "The rotateAxisLabels() method gets at most two input numeric values representing the X-axis and Y-axis label orientations from the horizontal axis, in degrees. " ...
                        + "If only one argument is provided, the input degree value will be used for both X and Y axis label orientations. " + newline ...
                        + "Alternatively, if no argument is passed, then all subplot axis labels will be rotated by the default 45 degrees." ...
                        + newline ...
                        );
            end

            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XLabel"); set(axisLabelX,"rotation",degreeX,"VerticalAlignment","top","HorizontalAlignment","right");
                    axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YLabel"); set(axisLabelY,"rotation",degreeY,"VerticalAlignment","middle","HorizontalAlignment","right");
                end
            end

        end % function rotateAxisLabels

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %***************************************************************************************************************************
    %***************************************************************************************************************************

    methods(Access=public, Hidden)

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function hideShow(self,varargin)

            allRequested = false;
            upperRequested = false;
            lowerRequested = false;
            methodName = varargin{end};
            methodNameIsShow = strcmp(methodName,"show");
            methodNameIsHide = strcmp(methodName,"hide");
            if methodNameIsShow
                switchValue = "on";
            elseif methodNameIsHide
                switchValue = "off";
            end
            if nargin==2
                allRequested = true;
            elseif nargin==3
                part = varargin{1};
                allRequested   = strcmpi(part,"all");
                diagRequested  = strcmpi(part,"diag");
                upperRequested = strcmpi(part,"upper");
                lowerRequested = strcmpi(part,"lower");
            else
                error   ( newline ...
                        + "The " + methodName + "() method gets only one input string value representing the ""part"" of the grid plot to " + methodName + "." + newline ...
                        + "Alternatively, if no argument is passed to " + methodName + "(), then all subplots will be affected." + newline ...
                        );
            end

            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    if allRequested || (icol<irow && lowerRequested) || (icol>irow && upperRequested) || (icol==irow && diagRequested)
                        set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"visible",switchValue);
                        set(get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"children"),"visible",switchValue);
                    end
                end
            end

            self.setAxisLabels();

        end % function hideShow

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function setAxisLabels(self)
        % show or hide axis labels and ticks depending on the presence of the neighbor subplots
            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    if irow < self.layout.axis.main.nrow
                        if strcmp(self.currentFig.subplotList{irow+1,icol}.currentFig.gca.Visible,"on")
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "XTick", [], "XTickLabel", []);
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.XLabel.String = "";
                        else
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "XTickMode", "auto", "XTickLabelMode", "auto");
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.XLabel.String = self.columns{icol};
                        end
                    end
                    if icol > 1
                        if strcmp(self.currentFig.subplotList{irow,icol-1}.currentFig.gca.Visible,"on")
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "YTick", [], "YTickLabel", [])
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLabel.String = "";
                        else
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "YTickMode", "auto", "YTickLabelMode", "auto")
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLabel.String = self.columns{irow};
                        end
                    end
                    degreeX = 45;
                    degreeY = 45;
                    if icol==1; degreeY = 90; end
                    if irow==self.layout.axis.main.nrow; degreeX = 0; end
                    axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XLabel"); set(axisLabelX,"rotation",degreeX,"VerticalAlignment","top","HorizontalAlignment","right");
                    axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YLabel"); set(axisLabelY,"rotation",degreeY,"VerticalAlignment","middle","HorizontalAlignment","right");
                end
            end
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function position = getColorbarAxisPosition(self)
            position =  [ self.layout.colorbar.origin.x ...
                        , self.layout.colorbar.origin.y ...
                        , self.layout.colorbar.width ...
                        , self.layout.colorbar.height ...
                        ];
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function position = getMainAxisPosition(self)
            position =  [ self.layout.axis.main.origin.x ...
                        , self.layout.axis.main.origin.y ...
                        , self.layout.axis.main.width ...
                        , self.layout.axis.main.height ...
                        ];
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function position = getSubplotAxisPosition(self,irow,icol)
            position =  [ (icol-1)*(self.layout.axis.subplot.interspace + self.layout.axis.subplot.width)   + self.layout.axis.main.margin.left ...
                        , (irow-1)*(self.layout.axis.subplot.interspace + self.layout.axis.subplot.height)  + self.layout.axis.main.margin.bottom ...
                        , self.layout.axis.subplot.width ...
                        , self.layout.axis.subplot.height ...
                        ];
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function updateLayout(self)
            self.layout.axis.main.height = 1 - self.layout.axis.main.margin.top;
            self.layout.axis.main.width = 1 - self.layout.axis.main.margin.right;
            self.layout.axis.main.nrow = length(self.columns);
            self.layout.axis.main.ncol = self.layout.axis.main.nrow;
            self.layout.axis.subplot.height = (1-self.layout.axis.main.margin.top-self.layout.axis.main.margin.bottom)/self.layout.axis.main.ncol - self.layout.axis.subplot.interspace;
            self.layout.axis.subplot.width  = (1-self.layout.axis.main.margin.left-self.layout.axis.main.margin.right)/self.layout.axis.main.nrow - self.layout.axis.subplot.interspace;
            self.layout.colorbar.origin.x = 1 - self.layout.axis.main.margin.right;
            self.layout.colorbar.origin.y = self.layout.axis.main.margin.bottom;
            self.layout.colorbar.height = self.layout.axis.main.nrow * ( self.layout.axis.subplot.height + self.layout.axis.subplot.interspace ) - self.layout.axis.subplot.interspace;
            self.layout.figWidth    = self.layout.figHeight *   ( 1 ...
                                                                - self.layout.axis.main.margin.bottom ...
                                                                + self.layout.axis.main.margin.right ...
                                                                + self.layout.axis.main.margin.left ...
                                                                - self.layout.axis.main.margin.top ...
                                                                + self.layout.colorbar.width ...
                                                                );
            if isfield(self,"currentFig") && isfield(self.currentFig,"gca") && isgraphics(self.currentFig.gca)
                set(self.currentFig.gca,"position",self.getMainAxisPosition());
            end
            if any(strcmp(fieldnames(self.currentFig),"colorbar"))
                set(self.currentFig.colorbar,"position",self.getColorbarAxisPosition());
            end
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef GridPlot