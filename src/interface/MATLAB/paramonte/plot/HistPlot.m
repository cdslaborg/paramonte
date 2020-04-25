classdef HistPlot < BasePlot

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        xcolumns
    end

    properties (Access = protected, Hidden)
        plotType
        isHist = false;
        isHist2 = false;
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
            if self.isHist
                self.histogram_kws = {};
            elseif self.isHist2
                self.histogram2_kws = {};
                self.colorbar_kws = {};
                self.ycolumns = {};
                self.colormap = {};
            elseif self.isHistfit
                self.histfit_kws = {};
            end
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

        function self = HistPlot(varargin) % expected input arguments: dataFrame, plotType
            self = self@BasePlot(varargin{1});
            self.plotType = lower(string(varargin{2}));
            if strcmp(self.plotType,"hist")
                self.isHist = true;
                prop="histogram_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif strcmp(self.plotType,"hist2")
                self.isHist2 = true;
                prop="histogram2_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="colorbar_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="ycolumns"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="colormap"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            elseif strcmp(self.plotType,"histfit")
                self.isHistfit = true;
                prop="histfit_kws"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
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

            self.parseArgs(varargin{:})

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            histEnabled = self.isHist && (isa(self.histogram_kws,"string") || isa(self.histogram_kws,"cell"));
            hist2Enabled = self.isHist2 && (isa(self.histogram2_kws,"string") || isa(self.histogram2_kws,"cell"));
            histfitEnabled = self.isHistfit && (isa(self.histfit_kws,"string") || isa(self.histfit_kws,"cell"));

            lgEnabled = ( isa(self.legend_kws,"string") || isa(self.legend_kws,"cell") );

            if histEnabled
                self.histogram_kws = addKeyVal("EdgeColor","none",self.histogram_kws{:});
            elseif hist2Enabled
                self.histogram2_kws = addKeyVal("EdgeColor","none",self.histogram2_kws{:});
                self.histogram2_kws = addKeyVal("FaceColor","flat",self.histogram2_kws{:});
                self.histogram2_kws = addKeyVal("DisplayStyle","bar3",self.histogram2_kws{:});
            %elseif histfitEnabled
                %self.histfit_kws = addKeyVal("EdgeColor","none",self.histfit_kws{:});
            end

            figEnabled = isa(self.gcf_kws,"string") || isa(self.gcf_kws,"cell");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % generate figure and axes if needed

            if figEnabled
                self.currentFig.gcf = figure(self.gcf_kws{:});
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

            if self.isHist2
                [ycolnames, ycolindex] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);
                ycolindexlen = length(ycolindex);
                maxLenColumns = max(maxLenColumns,ycolindexlen);
                if ycolindexlen~=maxLenColumns && ycolindexlen>1; error("length of ycolumns must be either 1 or equal to the maximum of the lengths of xcolumns."); end
                if ycolindexlen==1
                    ydata = self.dfref{rowindex,ycolindex};
                end
            end

            if xcolindexlen~=maxLenColumns && xcolindexlen>1; error("length of xcolumns must be either 1 or equal to the maximum of the lengths of ycolumns."); end

            % add line plot

            box on; grid on;

            lglabels = [];

            if histEnabled || hist2Enabled || histfitEnabled

                for i = 1:maxLenColumns

                    if xcolindexlen>1
                        xdata = self.dfref{rowindex,xcolindex(i)};
                    end
                    if hist2Enabled && ycolindexlen>1
                        ydata = self.dfref{rowindex,ycolindex(i)};
                    end

                    if lgEnabled
                        if hist2Enabled 
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

                    if histEnabled
                        self.currentFig.histogram = histogram   ( xdata ...
                                                                , self.histogram_kws{:} ...
                                                                );
                        hold on;
                    elseif hist2Enabled
                        self.currentFig.histogram2 = histogram2 ( xdata ...
                                                                , ydata ...
                                                                , self.histogram2_kws{:} ...
                                                                );
                        hold on;
                    elseif histfitEnabled
                        self.currentFig.histogram = histfit ( xdata ...
                                                            , self.histfit_kws{:} ...
                                                            );
                        hold on;
                    end

                end % loop plot

                self.currentFig.gca = gca;

            end

            if self.isHist2
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

            if self.isHist2
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

            colorbarEnabled = self.isHist2 && ( isa(self.colorbar_kws,"string") || isa(self.colorbar_kws,"cell") ); %&& ccolindexlen==1;
            if colorbarEnabled
                [fontsize, keyFound] = getKeyVal("fontsize",self.colorbar_kws{:});
                if keyFound
                    fontsize_kws = {"fontsize",fontsize};
                else
                    fontsize_kws = {"fontsize",self.currentFig.ylabel.FontSize};
                    self.colorbar_kws = [ self.colorbar_kws , fontsize_kws ];
                end
                self.currentFig.colorbar = colorbar(self.colorbar_kws{:});
                ylabel(self.currentFig.colorbar,"Density of Points",fontsize_kws{:});
            else
                colorbar('off');
                self.currentFig.colorbar = [];
            end

            if ~self.isHist2 || (self.isHist2 && ~isempty(self.legend_kws))
                self.doBasePlotStuff(lgEnabled,lglabels);
            end

            hold off;

        end % function plot

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef LinePlot