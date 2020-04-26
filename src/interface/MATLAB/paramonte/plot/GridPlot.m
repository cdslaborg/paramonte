classdef GridPlot < BasePlot

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = protected, Hidden)
        %dfref = [];
        %isdryrun = [];
        plotTypeList =  [ "line" ...
                        , "histogram" ...
                        , "histogram2" ...
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
        %histogram2_kws
        %histfit_kws
        %scatter_kws
        %scatter3_kws
        %lineScatter_kws
        %lineScatter3_kws
        colorbar_kws
        colormap
        layout
        target

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
                self.layout.subplot.(self.plotTypeList{i}) = struct();
            end
            self.colormap = {};
            
            self.colorbar_kws = struct();
            self.colorbar_kws.enabled = true;
            self.colorbar_kws.fontsize = [];

            self.layout.plotType.upper = "lineScatter";
            self.layout.plotType.lower = "histogram2";
            self.layout.plotType.diag = "histogram";
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

            self.isdryrun = true;
            self.plot();
            self.isdryrun = false;

            %self.target = GridTarget(self.dfref, self.columns, self.currentFig);

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

            self.gcf_kws.position = [ 50                    ... origin x
                                    , 100                   ... origin y
                                    , self.layout.figWidth  ... width
                                    , self.layout.figHeight ... height
                                    ];

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
                    self.layout.subplot.(self.plotTypeList(i)).rows = self.rows;
                else
                    self.layout.subplot.(self.plotTypeList(i)).rows = {};
                end

                %self.layout.subplot.(self.plotTypeList(i)).xcolumns = {};

                % set color properties

                if ~( strcmp(self.plotTypeList(i),"histogram") || strcmp(self.plotTypeList(i),"histfit") )

                    %self.layout.subplot.(self.plotTypeList(i)).ycolumns = {};

                    if ccolumnIsPresent
                        if contains(self.plotTypeList(i),"3")
                            self.layout.subplot.(self.plotTypeList(i)).zcolumns = self.ccolumn;
                        end
                        if strcmp(self.plotTypeList(i),"histogram2")
                            self.layout.subplot.(self.plotTypeList(i)).histogram2_kws.numbins = [100 100];
                            self.layout.subplot.(self.plotTypeList(i)).histogram2_kws.showemptybins = "off";
                            self.layout.subplot.(self.plotTypeList(i)).colormap = {};
                        else
                            self.layout.subplot.(self.plotTypeList(i)).ccolumns = self.ccolumn;
                            if ~isfield(self.layout.subplot.(self.plotTypeList(i)),"plot_kws") || isempty(self.layout.subplot.(self.plotTypeList(i)).plot_kws)
                                self.layout.subplot.(self.plotTypeList(i)).plot_kws = struct();
                            end
                            self.layout.subplot.(self.plotTypeList(i)).plot_kws.linewidth = 0.75;
                            self.layout.subplot.(self.plotTypeList(i)).plot_kws.color = uint8([200 200 200 100]);
                            self.layout.subplot.(self.plotTypeList(i)).plot_kws.enabled = true;
                            self.layout.subplot.(self.plotTypeList(i)).surface_kws.enabled = false;
                            self.layout.subplot.(self.plotTypeList(i)).colormap = {}; % self.colormap
                        end
                        self.layout.subplot.(self.plotTypeList(i)).colorbar_kws.enabled = false;
                    end

                end

                % check gcf_kws/box/grid presence

                self.layout.subplot.(self.plotTypeList(i)).legend_kws.enabled = false;
                self.layout.subplot.(self.plotTypeList(i)).gcf_kws.enabled = false;
                self.layout.subplot.(self.plotTypeList(i)).gca_kws.box = self.layout.axis.subplot.box;
                self.layout.subplot.(self.plotTypeList(i)).gca_kws.fontsize = self.layout.axis.subplot.fontsize;

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

            if self.gcf_kws.enabled
                gcf_kws_cell = convertStruct2Cell(self.gcf_kws,{"enabled","singleOptions"});
                if isfield(self.gcf_kws,"singleOptions"); gcf_kws_cell = { gcf_kws_cell{:} , self.gcf_kws.singleOptions{:} }; end
                self.currentFig.gcf = figure( gcf_kws_cell{:} );
            else
                self.currentFig.gcf = gcf;
            end
            set(0, "CurrentFigure", self.currentFig.gcf);

            % construct plot handles

            self.lowerEnabled = any(strcmp(self.plotTypeList,self.layout.plotType.lower));
            self.upperEnabled = any(strcmp(self.plotTypeList,self.layout.plotType.upper));
            self.diagEnabled = any(strcmp(self.plotTypeList,self.layout.plotType.diag));

            if self.lowerEnabled
                if isempty(self.layout.plotType.lower)
                    self.layout.plotType.lower = "histogram2";
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
                    self.layout.plotType.diag = "histogram";
                else
                    self.layout.plotType.diag = string(self.layout.plotType.diag);
                end
            end

            % generate axes and subplots

            self.currentFig.gca = axes  ( "position" ...
                                        , self.getMainAxisPosition()  ...
                                        ..., "Position", mainPosition ...
                                        , "Xlim", [0 1], "Ylim", [0 1] ...
                                        , "Color", "none" ...
                                        , "Parent", self.currentFig.gcf ...
                                        );
            axis(self.currentFig.gca,"off");

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
                diagPlotClass = "HistPlot";
            end

            % adjust limits

            columns.min = zeros(self.layout.axis.main.nrow,1);
            columns.max = zeros(self.layout.axis.main.nrow,1);
            for i = 1:self.layout.axis.main.nrow
                minvalue = min(self.dfref.(self.columns{i}));
                maxvalue = max(self.dfref.(self.columns{i}));
                margin = 0.1 * (maxvalue - minvalue);
                columns.min(i) = minvalue - margin;
                columns.max(i) = maxvalue + margin;
            end

            targetRGB = getRGB("deep carrot orange");

            iaxis = 0;
            self.currentFig.subplotList = cell(self.layout.axis.main.ncol,self.layout.axis.main.nrow);
            for irow = 1:self.layout.axis.main.nrow
                for icol = 1:self.layout.axis.main.ncol

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

                        xcolumns = self.columns(icol);
                        ycolumns = self.columns(irow);

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

                        plotType_kws = currentPlotType + "_kws";
                        fieldList = fieldnames(self.layout.subplot.(currentPlotType));
                        %valueList = struct2cell(self.layout.subplot.(currentPlotType));
                        for i = 1:length(fieldList)
                            if isa(self.layout.subplot.(currentPlotType).(fieldList{i}),"struct")
                                fieldSubList = fieldnames(self.layout.subplot.(currentPlotType).(fieldList{i}));
                                for j = 1:length(fieldSubList)
                                    self.currentFig.subplotList{irow,icol}.(fieldList{i}).(fieldSubList{j}) = self.layout.subplot.(currentPlotType).(fieldList{i}).(fieldSubList{j});
                                end
                            else
                                self.currentFig.subplotList{irow,icol}.(fieldList{i}) = self.layout.subplot.(currentPlotType).(fieldList{i});
                            end
                        end

                        if strcmp(currentPlotType,"histogram2")
                            %colormap();
                            %self.layout.subplot.(currentPlotType).colormap = flipud(cold());
                            %self.layout.subplot.(currentPlotType).colormap = "autumn";
                            %colormap(self.layout.subplot.(currentPlotType).colormap);
                        end

                        self.currentFig.subplotList{irow,icol}.plot()
                        self.currentFig.subplotList{irow,icol}.currentFig.gca = currentSubplotAxis;

                        for fname = ["hline_kws","vline_kws","scatter_kws"]
                            self.currentFig.subplotList{irow,icol}.target.(fname).color = targetRGB;
                        end
                        if strcmp(currentPlotType,"histogram")
                            self.currentFig.subplotList{irow,icol}.target.hline_kws.enabled = false;
                            self.currentFig.subplotList{irow,icol}.target.scatter_kws.enabled = false;
                        elseif strcmp(currentPlotType,"histogram2")
                            self.currentFig.subplotList{irow,icol}.target.hline_kws.enabled = false;
                            self.currentFig.subplotList{irow,icol}.target.vline_kws.enabled = false;
                            self.currentFig.subplotList{irow,icol}.target.scatter_kws.enabled = false;
                        end

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

                        self.currentFig.subplotList{irow,icol}.currentFig.gca.XLim = [ columns.min(icol) columns.max(icol) ];
                        if ~strcmp(currentPlotType,"histogram")
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLim = [ columns.min(irow) columns.max(irow) ];
                        end

                        %set(self.currentFig.axisList(iaxis)(iaxis),'box',self.layout.axis.subplot.box);
                        hold off;

                    else

                        disp( "skipping subplot #" + string(iaxis+": ("+irow+","+icol+")") + msgSuffix )
                        self.currentFig.subplotList{irow,icol} = [];

                    end

                    grid(self.currentFig.subplotList{irow,icol}.currentFig.gca,self.layout.axis.subplot.grid);

                end
            end

            % set up colorbar

            if self.colorbar_kws.enabled && cEnabled
                self.show("colorbar");
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

        end % function plot

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function hide(self,varargin)
            self.hideShow("hide",varargin{:})
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function show(self,varargin)
            self.hideShow("show",varargin{:})
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

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function addTarget(self,varargin)
        %   Add a target on an existing plot (the current active axes object)
        %   based on the 'value' attribute of the target object.
        %   
        %   Parameters
        %   ----------
        %       None
        %   
        %   Returns
        %   -------
        %       None. However, this method causes side-effects by manipulating 
        %       the existing attributes of the object.

            if nargin==1
                statistic = "mode";
            else
                statistic = [];
                for statistic = ["mode","mean","median"]
                    if strcmpi(varargin{1},statistic)
                        break;
                    end
                end
                if isempty(statistic)
                    [statistic, keyFound] = getKeyVal(key,varargin{:})
                    if ~keyFound
                        statistic = "mode";
                    end
                end
            end

            % parse args

            props = struct();
            for i = 1:length(varargin)
                for kws = ["hline_kws","vline_kws","scatter_kws","values"]
                    [val, keyFound] = getKeyVal(kws,varargin{:});
                    if keyFound
                        props.(kws) = val;
                    end
                end
            end

            % compute the statistics

            columnsLen = length(self.columns);
            if isfield(props,"values")
                values = props.values;
            else
                values = zeros(columnsLen,1);
                statIsMode = strcmpi(statistic,"mode");
                statIsMean = strcmpi(statistic,"mean");
                statIsMedian = strcmpi(statistic,"median");
                if statIsMode
                    [~,irowMax] = max(self.dfref.SampleLogFunc);
                end
                for i = 1:columnsLen
                    if statIsMode
                        values(i) = self.dfref.(self.columns{i})(irowMax);
                    elseif statIsMean
                        values(i) = mean(self.dfref.(self.columns{i}));
                    elseif statIsMedian
                        values(i) = median(self.dfref.(self.columns{i}));
                    end
                end
            end

            for icol = 1:columnsLen
                for irow = 1:columnsLen

                    set(self.currentFig.gcf, "CurrentAxes", self.currentFig.subplotList{irow,icol}.currentFig.gca);

                    hold on;

                    self.currentFig.subplotList{irow,icol}.target.value = [ values(icol) values(irow) ];
                    for kws = ["hline_kws","vline_kws","scatter_kws"]
                        if isfield(props,kws)
                            for i = 1:2:length(props.(kws))
                                key = props.(kws){i};
                                val = props.(kws){i+1};
                                self.currentFig.subplotList{irow,icol}.target.(kws).(key) = val;
                            end
                        end
                    end
                    self.currentFig.subplotList{irow,icol}.target.plot();

                    hold off;

                end
            end

        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %***************************************************************************************************************************
    %***************************************************************************************************************************

    methods(Access=public, Hidden)

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function hideShow(self,varargin)

            diagRequested = false;
            lowerRequested = false;
            upperRequested = false;
            colorbarRequested = false;
            methodName = varargin{1};
            methodNameIsShow = strcmp(methodName,"show");
            methodNameIsHide = strcmp(methodName,"hide");
            if methodNameIsShow
                switchValue = "on";
            elseif methodNameIsHide
                switchValue = "off";
            end
            if nargin==2 % act on all parts of the figure
                diagRequested = true;
                lowerRequested = true;
                upperRequested = true;
                colorbarRequested = true;
            elseif nargin>2
                for i = 2:length(varargin)
                    part = varargin{i};
                    if strcmpi(part,"diag")
                        diagRequested  = true;
                    elseif strcmpi(part,"upper")
                        upperRequested = true;
                    elseif strcmpi(part,"lower")
                        lowerRequested = true;
                    elseif strcmpi(part,"colorbar")
                        colorbarRequested = true;
                    else
                        error   ( newline ...
                                + "Unrecognized input argument pass to the " + methodName + "() method: " + part ...
                                + "The " + methodName + "() method accepts only a list of comma-separated string values representing the ""part"" of the grid plot to " + methodName + ". " ...
                                + "Alternatively, if no argument is passed to " + methodName + "(), then all subplots will be affected." + newline ...
                                );
                    end
                end
            else
                error   ( newline ...
                        + "Internal error occurred. Not enough input arguments to hideShow()" ...
                        + newline ...
                        );
            end

            if colorbarRequested

                if methodNameIsShow

                    colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});

                    colormap(self.colormap);

                    ccolumnValues = self.dfref.(self.ccolumn);
                    colorbarLimit = [ min(ccolumnValues) max(ccolumnValues) ]; % Colorbar range
                    caxis(colorbarLimit);
                    set(self.currentFig.gca,"CLim",colorbarLimit);
                    self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});

                    colorbarLabel = self.ccolumn;
                    ylabel(self.currentFig.colorbar,colorbarLabel,"fontsize",self.layout.colorbar.fontsize);

                    set ( self.currentFig.colorbar ...
                        , "position" ...
                        , self.getColorbarAxisPosition() ...
                        , "FontSize", self.layout.colorbar.fontsize ...
                        ) ;

                elseif methodNameIsHide

                    colorbar("off");

                end

            end

            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    if (icol<irow && lowerRequested) || (icol>irow && upperRequested) || (icol==irow && diagRequested)
                        set(self.currentFig.subplotList{irow,icol}.currentFig.gca,"visible",switchValue);
                        set(get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"children"),"visible",switchValue);
                    end
                end
            end

            if lowerRequested || upperRequested || diagRequested
                self.setAxisLabels();
            end

        end % function hideShow

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function setAxisLabels(self)
        % show or hide axis labels and ticks depending on the presence of the neighbor subplots
            for icol = 1:self.layout.axis.main.ncol
                for irow = 1:self.layout.axis.main.nrow
                    if irow < self.layout.axis.main.nrow
                        if strcmp(self.currentFig.subplotList{irow+1,icol}.currentFig.gca.Visible,"on")
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "XTickLabel", []);
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.XLabel.String = "";
                        else
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "XTickLabelMode", "auto");
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.XLabel.String = self.columns{icol};
                        end
                    end
                    if icol > 1
                        if strcmp(self.currentFig.subplotList{irow,icol-1}.currentFig.gca.Visible,"on")
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "YTickLabel", [])
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLabel.String = "";
                        else
                            set(self.currentFig.subplotList{irow,icol}.currentFig.gca, "YTickLabelMode", "auto")
                            self.currentFig.subplotList{irow,icol}.currentFig.gca.YLabel.String = self.columns{irow};
                        end
                    end
                    if ~(icol==1)
                        axisLabelY = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"YLabel"); set(axisLabelY,"rotation",45,"VerticalAlignment","middle","HorizontalAlignment","right");
                    end
                    if ~(irow==self.layout.axis.main.nrow)
                        axisLabelX = get(self.currentFig.subplotList{irow,icol}.currentFig.gca,"XLabel"); set(axisLabelX,"rotation",45,"VerticalAlignment","top","HorizontalAlignment","right");
                    end
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
            if isfield(self.currentFig,"subplotList")
                for irow = 1:self.layout.axis.main.nrow
                    for icol = 1:self.layout.axis.main.ncol
                        jrow = self.layout.axis.main.nrow-irow+1;
                        set ( self.currentFig.subplotList{irow,icol}.currentFig.gca, "position", self.getSubplotAxisPosition(jrow,icol) );
                    end
                end
            end
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef GridPlot