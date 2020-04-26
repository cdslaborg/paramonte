classdef LineScatterPlot < BasePlot
%   This is the LinePlot class for generating instances 
%   of line figures based on matplotlib library's 
%   line() and functions.
%
%   Usage
%   -----
%   first generate an object of this class by optionally 
%   passing the following parameters described below. Then call 
%   the plot() method. The generated object is also callable with 
%   the same input parameters as the object's constructor.
%
%   Parameters
%   ----------
%       dataFrame
%           a Pandas dataframe from which the selected data will be plotted.
%       xcolumns
%           optional argument that determines the columns of dataframe to serve as 
%           the x-values. It can have three forms:
%               1.  a list of column indices in the input dataFrame.
%               2.  a list of column names in dataFrame.columns.
%               3.  a range(start,stop,step), representing the column indices 
%                   in the input dataFrame.
%           Examples:
%               1.  xcolumns = [0,1,4,3]
%               2.  xcolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  xcolumns = range(17,7,-2)
%           However, in all cases, it must have a length that is either 1 or equal 
%           to the length of ycolumns. If the length is 1, then xcolumns will be 
%           plotted against data corresponding to each element of ycolumns.
%           If not provided, the default will be the count of the rows of the 
%           input dataFrame.
%       ycolumns
%           optional argument that determines the columns of dataframe to serve 
%           as the y-values. It can have three forms:
%               1.  a list of column indices in the input dataFrame.
%               2.  a list of column names in dataFrame.columns.
%               3.  a range(start,stop,step), representing the column indices 
%                   in the input dataFrame.
%           Examples:
%               1.  ycolumns = [0,1,4,3]
%               2.  ycolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  ycolumns = range(17,7,-2)
%           However, in all cases, it must have a length that is either 1 or 
%           equal to the length of xcolumns. If the length is 1, then ycolumns 
%           will be plotted against data corresponding to each element of xcolumns. 
%           If not provided, the default includes all columns of the dataframe. 
%       ccolumns
%           (line-color-columns) optional argument that determines the columns 
%           of dataframe to serve as the color-values corresponding to each 
%           line-segment in the plot. If provided, a matplotlib LineCollection 
%           will be added to the plot. It can have three forms:
%               1.  a list of column indices in the input dataFrame.
%               2.  a list of column names in dataFrame.columns.
%               3.  a range(start,stop,step), representing the column indices 
%                   in the input dataFrame.
%           Examples:
%               1.  ccolumns = [0,1,4,3]
%               2.  ccolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  ccolumns = range(17,7,-2)
%           However, in all cases, it must have a length that is either 1 or 
%           equal to the lengths of xcolumns or ycolumns, whichever is not 1. 
%           If the length is 1, then the same color will be used for plotting 
%           data corresponding to each element of xcolumns. If not provided, 
%           the default will be the count of the rows of the input dataFrame. 
%           If set to None, no colored LineCollection will be added to the plot.
%       rows
%           optional argument that determines the rows of dataframe 
%           to be visualized. It can be either:
%               1.  a range(start,stop,step), or, 
%               2.  a list of row indices in dataFrame.index.
%           Examples:
%               1.  rows = range(17,7,-2)
%               2.  rows = [i for i in range(7,17)]
%           If not provided, the default includes all rows of the dataframe.
%       lc_kws
%           optional cell array of string values to be passed to matplotlib's 
%           LineCollection() class. For example: 
%               lc_kws = {"cmap":"autumn"}
%           The default is {}. 
%           If set to None, no colored LineCollection will be plotted.
%       gca_kws
%           optional cell array of string values to be passed to seaborn's 
%           set() function. For example: 
%               gca_kws = {"style":"darkgrid"}
%           if set to None, then no call to set() function will be made.
%       gcf_kws
%           optional cell array of string values to be passed to matplotlib's 
%           figure() function. For example: 
%               gcf_kws = {"facecolor":"w","dpi":150}
%       plot_kws
%           optional cell array of string values to be passed to matplotlib's 
%           plot() function. For example: 
%               plot_kws = {"linewidth":0.5}
%           The default is {}.
%           If set to None, no fixed-color line plot will be plotted. 
%       legend_kws
%           optional cell array of string values to be passed to matplotlib's 
%           legend() function. If it is set to None, no legend will be added
%           to the plot. If it is set to {}, then the legend will be added 
%           to the plot and automatically adjusted. For example: 
%               legend_kws = {"labels":["Variable1","Variable2"]}
%           legend will be added to plot only if simple line plot with no 
%           color-mappings are requested.
%       colorbar_kws
%           optional cell array of string values to be passed to matplotlib's 
%           figure.colorbar() function. For example: 
%               colorbar_kws = {"orientation":"vertical"}
%           The default is {}. If set to None, no colorbar will be plotted.
%       outputFile
%           optional string representing the name of the output file in which 
%           the figure will be saved. If not provided, no output file will be generated.
%       axes
%           the axes object to which the plot must be added.
%           The default is None in which case the output of matplotlib's gca()
%           function wil be used to get the current active axes.
%
%   Attributes
%   ----------
%       All of the parameters described above, except dataFrame.
%           a reference to the dataFrame will be implicitly stored in the object.
%       target
%           a callable object of ParaMonte library's class Target which can 
%           be used to add target point and lines to the current active plot.
%       currentFig
%           an object of class ParaMonteFigure which is initially None, but upon 
%           making a plot, is populated with attributes representing all outputs 
%           from matplotlib/seaborn function calls, with the following attributes:
%               figure
%                   the output of matplotlib's figure() function.
%               axes
%                   the axes on which the plot is made.
%               plot
%                   a list of Line2D objects which is the output 
%                   of matplotlib's plot() function.
%               lineCollection
%                   an object of type LineCollection which is the 
%                   output of matplotlib's LineCollection() class.
%               legend
%                   the output of matplotlib's legend() function.
%               colorbar
%                   the output of matplotlib's Figure.colorbar() function.
%
%   Returns
%   -------
%       LinePlot
%
%   ----------------------------------------------------------------------

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)

        ccolumns
        colorbar_kws
        colormap
        target

    end

    properties (Access = protected, Hidden)
        %dfref = [];
        %isdryrun = [];
        is3d
        plotType
        isLinePlot = false;
        isScatterPlot = false;
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function reset(self)

            reset@BasePlot(self);
            prop="xcolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            self.xcolumns = {};
            self.ycolumns = {};
            self.ccolumns = {};
            self.colormap = {};

            self.isLinePlot = false;
            if contains(self.plotType,"line")
                self.isLinePlot = true;
                prop="plot_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="surface_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end
            self.isScatterPlot = false;
            if contains(self.plotType,"scatter")
                self.isScatterPlot = true;
                prop="scatter_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end
            if contains(self.plotType,"3")
                prop="zcolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.is3d = true;
            else
                self.is3d = false;
            end

            if self.isLinePlot

                self.plot_kws = struct();
                self.plot_kws.enabled = true;
                self.plot_kws.linewidth = {};
                self.plot_kws.singleOptions = {};
                self.plot_kws.color = {};

                self.surface_kws = struct();
                self.surface_kws.enabled = true;
                self.surface_kws.facecolor = "none";
                self.surface_kws.edgecolor = "flat";
                self.surface_kws.edgealpha = 0.5;
                self.surface_kws.linestyle = "-";
                self.surface_kws.marker = "none";

            end

            if self.isScatterPlot
                self.scatter_kws = struct();
                self.scatter_kws.marker = [];
                self.scatter_kws.singleOptions = {};
                self.scatter_kws.size = [];
                self.scatter_kws.enabled = true;
            end

            self.colorbar_kws = struct();
            self.colorbar_kws.enabled = true;
            self.colorbar_kws.fontsize = [];

            self.isdryrun = true;
            self.plot();
            self.isdryrun = false;

            self.target = Target();

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = LineScatterPlot(varargin) % expected input arguments: dataFrame, plotType

            self = self@BasePlot(varargin{1});
            self.plotType = lower(string(varargin{2}));
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

            self.parseArgs(varargin{:})

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %if self.scatter_kws.enabled
            %    if getKeyVal(self.scatter_kws)
            %    key = string(self.scatter_kws)
            %    if strcmpcolorself.scatter_kws.enabled = 
            %end

            cEnabled =  ( isa(self.colormap,"string") || isa(self.colormap,"cell") || isa(self.colormap,"char") ) && ...
                        ( isa(self.ccolumns,"string") || isa(self.ccolumns,"cell") || isa(self.ccolumns,"char") ) && ...
                        ( ( self.isLinePlot && self.surface_kws.enabled ) || ( self.isScatterPlot && self.scatter_kws.enabled ) );
            lgEnabled = self.legend_kws.enabled && ~cEnabled;

            %if any(strcmp(properties(self),"scatter_kws")) && self.scatter_kws.enabled
            if self.isScatterPlot && self.scatter_kws.enabled
                key = "marker"; val = ".";
                if isfield(self.scatter_kws,key) && isempty(self.scatter_kws.(key))
                    self.scatter_kws.(key) = val;
                end
                key = "size"; val = 10;
                if isfield(self.scatter_kws,key) && isempty(self.scatter_kws.(key))
                    self.scatter_kws.(key) = val;
                end
            end
%if self.isLinePlot
%"begin self.plot_kws"
%self.plot_kws
%"end self.plot_kws"
%end
%return
            if self.isLinePlot && self.plot_kws.enabled
                key = "linewidth"; val = 1;
                if isfield(self.plot_kws,key) && isempty(self.plot_kws.(key))
                    self.plot_kws.(key) = val;
                end
            end

            if cEnabled && ~getVecLen(self.colormap)
                if self.is3d
                    self.colormap = "winter";
                else
                    self.colormap = "autumn";
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
                set(0, "CurrentFigure", gcf);
                self.currentFig.gcf = gcf;
            end

            % check rows presence

            if getVecLen(self.rows)
                rowindex = self.rows;
            else
                rowindex = 1:1:length(self.dfref{:,1});
            end
            %ndata = length(rowindex);

            % check columns presence

            if getVecLen(self.xcolumns)
                [xcolnames, xcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.xcolumns);
            else
                xcolindex = [];
                xcolnames = "Count";
                xdata = rowindex;
            end

            [ycolnames, ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);

            if self.is3d && getVecLen(self.zcolumns)
                [zcolnames, zcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.zcolumns);
            else
                zcolindex = [];
                zcolnames = "Count";
                zdata = rowindex;
            end

            % set color data

            if cEnabled
                if getVecLen(self.ccolumns)
                    [ccolnames, ccolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ccolumns);
                else
                    ccolindex = [];
                    ccolnames = "Count";
                    cdata = 1:1:length(self.dfref{:,1});
                end
            else
                ccolindex = [];
                ccolnames = [];
                cdata = [];
            end

            % check the lengths are consistent

            xcolindexlen = length(xcolindex);
            ycolindexlen = length(ycolindex);
            zcolindexlen = length(zcolindex);
            ccolindexlen = length(ccolindex);
            maxLenColumns = max (   [ xcolindexlen ...
                                    , ycolindexlen ...
                                    , zcolindexlen ...
                                    , ccolindexlen ...
                                    ] ...
                                );

            if xcolindexlen~=maxLenColumns && xcolindexlen>1; error("length of xcolumns must be either 1 or equal to the maximum of the lengths of ycolumns, zcolumns, ccolumns."); end
            if ycolindexlen~=maxLenColumns && ycolindexlen>1; error("length of ycolumns must be either 1 or equal to the maximum of the lengths of xcolumns, zcolumns, ccolumns."); end
            if zcolindexlen~=maxLenColumns && zcolindexlen>1; error("length of zcolumns must be either 1 or equal to the maximum of the lengths of xcolumns, ycolumns, ccolumns."); end
            if ccolindexlen~=maxLenColumns && ccolindexlen>1; error("length of ccolumns must be either 1 or equal to the maximum of the lengths of xcolumns, ycolumns, zcolumns."); end

            % assign data in case of single column assignments

            if xcolindexlen==1
                xdata = self.dfref{rowindex,xcolindex};
            end
            if ycolindexlen==1
                ydata = self.dfref{rowindex,ycolindex};
            end
            if zcolindexlen==1
                zdata = self.dfref{rowindex,zcolindex};
            end
            if ccolindexlen==1
                cdata = self.dfref{rowindex,ccolindex};
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%

            if self.isLinePlot
                plot_kws_cell = convertStruct2Cell(self.plot_kws,{"enabled","singleOptions"});
            end
            if self.isScatterPlot
                scatter_kws_cell = convertStruct2Cell(self.scatter_kws,{"enabled","singleOptions","color","size"});
            end
            if self.isLinePlot
                surface_kws_cell = convertStruct2Cell(self.surface_kws,{"enabled","singleOptions","color","size"});
            end

            %%%%%%%%%%%%%%%
            % add line plot
            %%%%%%%%%%%%%%%

            box on; grid on;

            lglabels = [];
            if cEnabled
                colormap(self.colormap);
            end

            if (self.isScatterPlot && self.scatter_kws.enabled) || (self.isLinePlot && (self.surface_kws.enabled || self.plot_kws.enabled))

                if ~self.is3d
                    zdata = zeros(length(rowindex(:)),1);
                end

                for i = 1:maxLenColumns

                    if xcolindexlen>1
                        xdata = self.dfref{rowindex,xcolindex(i)};
                    end
                    if ycolindexlen>1
                        ydata = self.dfref{rowindex,ycolindex(i)};
                    end
                    if zcolindexlen>1
                        zdata = self.dfref{rowindex,zcolindex(i)};
                    end
                    if ccolindexlen>1
                        cdata = self.dfref{rowindex,ccolindex(i)};
                    end

                    if lgEnabled && ~self.is3d
                        if xcolindexlen<2 && ycolindexlen>=1
                            lglabels = [ lglabels , ycolnames(i) ];
                        elseif xcolindexlen>1 && ycolindexlen<2
                            lglabels = [ lglabels , xcolnames(i) ];
                        else
                            lglabels = [ lglabels , xcolnames(i)+"-"+ycolnames(i) ];
                        end
                    end

                    % add plot

                    if self.isLinePlot
                        if self.plot_kws.enabled
                            if self.is3d
                                if self.surface_kws.enabled
                                    self.currentFig.surface = surface   ( "XData",[xdata(:) xdata(:)] ...
                                                                        , "YData",[ydata(:) ydata(:)] ...
                                                                        , "ZData",[zdata(:) zdata(:)] ...
                                                                        , "CData",[cdata(:) cdata(:)] ...
                                                                        , surface_kws_cell{:} ...
                                                                        );
                                    view(3);
                                    hold on;
                                else
                                    self.currentFig.plot3 = plot3   ( xdata ...
                                                                    , ydata ...
                                                                    , zdata ...
                                                                    , plot_kws_cell{:} ...
                                                                    );
                                end
                            else
                                self.currentFig.plot = plot ( xdata ...
                                                            , ydata ...
                                                            , plot_kws_cell{:} ...
                                                            );
                            end
                            hold on;
                        end
                    end

                    if self.isScatterPlot && self.scatter_kws.enabled
                        if cEnabled
                            if self.is3d
                                self.currentFig.scatter3 = scatter3 ( xdata ...
                                                                    , ydata ...
                                                                    , zdata ...
                                                                    , self.scatter_kws.size ...
                                                                    , cdata ...
                                                                    , scatter_kws_cell{:} ...
                                                                    );
                            else
                                self.currentFig.scatter = scatter   ( xdata ...
                                                                    , ydata ...
                                                                    , self.scatter_kws.size ...
                                                                    , cdata ...
                                                                    , scatter_kws_cell{:} ...
                                                                    );
                            end
                        else
                            if self.is3d
                                %plot_kws = {};
                                %if ~isa(self.plot_kws,"cell"); plot_kws = self.plot_kws;
                                self.currentFig.plot = plot3( xdata ...
                                                            , ydata ...
                                                            , zdata ...
                                                            , scatter_kws_cell{:} ...
                                                            );
                            else
                                %plot_kws = {};
                                %if ~isa(self.plot_kws,"cell"); plot_kws = self.plot_kws;
                                self.currentFig.plot = plot ( xdata ...
                                                            , ydata ...
                                                            , scatter_kws_cell{:} ...
                                                            );
                            end
                        end
                        hold on;
                    end

                end % loop plot

                self.currentFig.gca = gca;

            end

            % add axis labels

            if xcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values");
            else
                self.currentFig.xlabel = xlabel(xcolnames(1));
            end

            if ycolindexlen>1
                self.currentFig.ylabel = ylabel("Variable Values");
            else
                self.currentFig.ylabel = ylabel(ycolnames(1));
            end

            if self.is3d
            if zcolindexlen>1
                self.currentFig.zlabel = zlabel("Variable Values");
            else
                self.currentFig.zlabel = zlabel(zcolnames(1));
            end
            end

            % add line colorbar

            if self.colorbar_kws.enabled && cEnabled
                if isempty(self.colorbar_kws.fontsize) || ~isa(self.colorbar_kws.fontsize,"numeric")
                    self.colorbar_kws.fontsize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,ccolnames(1),"fontsize",self.colorbar_kws.fontsize);
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if ~self.is3d || (self.is3d && self.legend_kws.enabled)
                self.doBasePlotStuff(lgEnabled,lglabels)
            end

            grid on; hold off;

        end % function plot

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef LinePlot