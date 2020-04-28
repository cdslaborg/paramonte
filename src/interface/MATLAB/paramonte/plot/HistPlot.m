classdef HistPlot < BasePlot

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        xcolumns
        target
    end

    properties (Access = protected, Hidden)
        plotType
        isHistogram = false;
        isHistogram2 = false;
        isHistfit = false;
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function reset(self)

            reset@BasePlot(self);
            self.xcolumns = {};

            if self.isHistogram

                self.histogram_kws = struct();
                self.histogram_kws.enabled = true;
                self.histogram_kws.edgecolor = {};
                self.histogram_kws.singleOptions = {};

            elseif self.isHistogram2

                self.ycolumns = {};
                self.colormap = {};

                self.histogram2_kws = struct();
                self.histogram2_kws.enabled = true;
                self.histogram2_kws.edgecolor = {};
                self.histogram2_kws.facecolor = {};
                self.histogram2_kws.displaystyle = {};
                self.histogram2_kws.singleOptions = {};

                self.colorbar_kws = struct();
                self.colorbar_kws.enabled = true;
                self.colorbar_kws.fontsize = [];
                self.colorbar_kws.singleOptions = {};

                self.legend_kws.enabled = false;

            elseif self.isHistfit

                self.histfit_kws = struct();
                self.histfit_kws.enabled = true;
                self.histfit_kws.singleOptions = {};

            end

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

        function self = HistPlot(varargin) % expected input arguments: dataFrame, plotType
            self = self@BasePlot(varargin{1});
            self.plotType = lower(string(varargin{2}));
            if strcmp(self.plotType,"histogram")
                self.isHistogram = true;
                prop="histogram_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif strcmp(self.plotType,"histogram2")
                self.isHistogram2 = true;
                prop="histogram2_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="colorbar_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="colormap"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif strcmp(self.plotType,"histfit")
                self.isHistfit = true;
                prop="histfit_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end
            self.reset();
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
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.isHistogram && self.histogram_kws.enabled
                fname = "histogram_kws";
                key = "edgecolor"; val = "none"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
            end

            if self.isHistogram2 && self.histogram2_kws.enabled
                fname = "histogram2_kws";
                key = "edgecolor"; val = "none"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "facecolor"; val = "flat"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
                key = "displaystyle"; val = "bar3"; if isfield(self.(fname),key) && isempty(self.(fname).(key)); self.(fname).(key) = val; end
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

            % check columns presence

            if getVecLen(self.xcolumns)
                [xcolnames, xcolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.xcolumns);
            else
                xcolindex = [];
                xcolnames = "Count";
                xdata = rowindex;
            end
            xcolindexlen = length(xcolindex);
            maxLenColumns = xcolindexlen;
            if xcolindexlen==1
                xdata = self.dfref{rowindex,xcolindex};
            end

            if self.isHistogram2
                [ycolnames, ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);
                ycolindexlen = length(ycolindex);
                maxLenColumns = max(maxLenColumns,ycolindexlen);
                if ycolindexlen~=maxLenColumns && ycolindexlen>1; error("length of ycolumns must be either 1 or equal to the maximum of the lengths of xcolumns."); end
                if ycolindexlen==1
                    ydata = self.dfref{rowindex,ycolindex};
                end
            end

            if xcolindexlen~=maxLenColumns && xcolindexlen>1; error("length of xcolumns must be either 1 or equal to the maximum of the lengths of ycolumns."); end

            %%%%%%%%%%%%%%%%%%%%%%%
            % get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%

            if self.isHistogram
                histogram_kws_cell = convertStruct2Cell(self.histogram_kws,{"enabled","singleOptions"});
            end
            if self.isHistogram2
                histogram2_kws_cell = convertStruct2Cell(self.histogram2_kws,{"enabled","singleOptions"});
            end
            if self.isHistfit
                histfit_kws_cell = convertStruct2Cell(self.histfit_kws,{"enabled","singleOptions"});
            end

            %%%%%%%%%%%%%%%
            % add line plot
            %%%%%%%%%%%%%%%

            box on; grid on;

            lglabels = [];

            if (self.isHistogram && self.histogram_kws.enabled) || (self.isHistogram2 && self.histogram2_kws.enabled) || (self.isHistfit && self.histfit_kws.enabled)

                for i = 1:maxLenColumns

                    if xcolindexlen>1
                        xdata = self.dfref{rowindex,xcolindex(i)};
                    end
                    if (self.isHistogram2 && self.histogram2_kws.enabled) && ycolindexlen>1
                        ydata = self.dfref{rowindex,ycolindex(i)};
                    end

                    if self.legend_kws.enabled
                        if (self.isHistogram2 && self.histogram2_kws.enabled) 
                            if xcolindexlen<2 && ycolindexlen>=1
                                lglabels = [ lglabels , ycolnames(i) ];
                            elseif xcolindexlen>1 && ycolindexlen<2
                                lglabels = [ lglabels , xcolnames(i) ];
                            else
                                lglabels = [ lglabels , xcolnames(i)+"-"+ycolnames(i) ];
                            end
                        else
                            lglabels = [ lglabels , xcolnames(i) ];
                        end
                    end

                    % add plot

                    if self.isHistogram && self.histogram_kws.enabled
                        self.currentFig.histogram = histogram   ( xdata ...
                                                                , histogram_kws_cell{:} ...
                                                                );
                        hold on;
                    elseif self.isHistogram2 && self.histogram2_kws.enabled
                        self.currentFig.histogram2 = histogram2 ( xdata ...
                                                                , ydata ...
                                                                , histogram2_kws_cell{:} ...
                                                                );
                        hold on;
                    elseif self.isHistfit && self.histfit_kws.enabled
                        self.currentFig.histogram = histfit ( xdata ...
                                                            , histfit_kws_cell{:} ...
                                                            );
                        hold on;
                    end

                end % loop plot

                self.currentFig.gca = gca;

            end

            if self.isHistogram2
                if ~getVecLen(self.colormap)
                    self.colormap = flipud(cold());
                end
                colormap(self.colormap);
            end

            % add axis labels

            if xcolindexlen>1
                self.currentFig.xlabel = xlabel("Variable Values");
            else
                self.currentFig.xlabel = xlabel(xcolnames(1));
            end

            if self.isHistogram2
                if ycolindexlen>1
                    self.currentFig.ylabel = ylabel("Variable Values");
                else
                    self.currentFig.ylabel = ylabel(ycolnames(1));
                end
                self.currentFig.zlabel = zlabel("Count");
            else
                self.currentFig.ylabel = ylabel("Count");
            end

            % add line colorbar

            colorbarEnabled = self.isHistogram2 && self.colorbar_kws.enabled;
            if colorbarEnabled
                if isempty(self.colorbar_kws.fontsize) || ~isa(self.colorbar_kws.fontsize,"numeric")
                    self.colorbar_kws.fontsize = self.currentFig.ylabel.FontSize;
                end
                colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","singleOptions"});
                self.currentFig.colorbar = colorbar(colorbar_kws_cell{:});
                ylabel(self.currentFig.colorbar,"Density of Points","fontsize",self.colorbar_kws.fontsize);
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            self.doBasePlotStuff(self.legend_kws.enabled,lglabels);

            hold off;

        end % function plot

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef HistPlot