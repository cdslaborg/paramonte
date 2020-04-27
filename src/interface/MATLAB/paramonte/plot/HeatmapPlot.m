classdef HeatmapPlot < BasePlot

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        title
        xcolumns
        ycolumns
        colormap
        precision
        heatmap_kws
        colorbar_kws
    end

    properties (Access = protected, Hidden)
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function reset(self)

            reset@BasePlot(self);
            self.xcolumns = {};
            self.ycolumns = {};

            self.precision = [];

            title = [];
            self.heatmap_kws = struct();
            self.heatmap_kws.enabled = true;
            self.heatmap_kws.ColorLimits = [];
            self.heatmap_kws.singleOptions = {};

            self.colorbar_kws = struct();
            %self.colorbar_kws.label = [];
            self.colorbar_kws.enabled = true;
            %self.colorbar_kws.fontsize = [];
            %self.colorbar_kws.singleOptions = {};

            self.colormap = [];

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

        function self = HeatmapPlot(dataFrame) % expected input arguments: dataFrame
            self = self@BasePlot(dataFrame,"heatmap");
            self.reset();
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function adjustColorLimits(self,newLimits)
            if nargin==1
                dum = max(abs(self.currentFig.heatmap.ColorLimits));
                newLimits = [-dum dum];
            elseif nargin~=2
                error   ( newline ...
                        + "colorLimits() method takes only one optional argument newLimits. If newLimits is provided as input argument, "  ...
                        + "the ColorLimits property of the heatmap plot will be set to newLimits." ...
                        + newline ...
                        );
            end
            self.currentFig.heatmap.ColorLimits = newLimits;
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

            parseArgs(self,varargin{:})


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if getVecLen(self.xcolumns)
                [xcolnames, ~] = getColNamesIndex(self.dfref.Properties.VariableNames,self.xcolumns);
            else
                self.dfref
                self.xcolumns = self.dfref.Properties.VariableNames;
            end
            if getVecLen(self.ycolumns)
                [ycolnames, ~] = getColNamesIndex(self.dfref.Properties.VariableNames,self.ycolumns);
            else
                self.ycolumns = self.dfref.Properties.RowNames;
            end

            % generate figure and axes if needed

            if self.gcf_kws.enabled
                gcf_kws_cell = convertStruct2Cell(self.gcf_kws,{"enabled","singleOptions"});
                if isfield(self.gcf_kws,"singleOptions"); gcf_kws_cell = { gcf_kws_cell{:} , self.gcf_kws.singleOptions{:} }; end
                self.currentFig.gcf = figure( gcf_kws_cell{:} );
            else
                set(0, "CurrentFigure", gcf);
                self.currentFig.gcf = gcf;
            end

            %%%%%%%%%%%%%%%%%%%%%%%
            % get keyword arguments
            %%%%%%%%%%%%%%%%%%%%%%%

            heatmap_kws_cell = convertStruct2Cell(self.heatmap_kws,{"enabled","ColorLimits","singleOptions"});

            %%%%%%%%%%%%%
            % add heatmap
            %%%%%%%%%%%%%

            if self.heatmap_kws.enabled
                if isempty(self.precision)
                    self.currentFig.heatmap = heatmap( xcolnames, ycolnames, self.dfref{xcolnames,ycolnames} );
                elseif isa(self.precision,"numeric")
                    self.currentFig.heatmap = heatmap( xcolnames, ycolnames, round(self.dfref{xcolnames,ycolnames},self.precision) );
                end
            end

            self.currentFig.gca = gca;

            if ~isempty(self.heatmap_kws.ColorLimits)
                self.currentFig.heatmap.ColorLimits = self.heatmap_kws.ColorLimits;
            end

            % add line colorbar

            if ~isempty(self.colormap)
                colormap(self.colormap);
            else
                colormap(redblue());
            end

            if ~self.colorbar_kws.enabled
                %if isempty(self.colorbar_kws.fontsize) || ~isa(self.colorbar_kws.fontsize,"numeric")
                %    self.colorbar_kws.fontsize = self.currentFig.gca.FontSize;
                %end
                %colorbar_kws_cell = convertStruct2Cell(self.colorbar_kws,{"enabled","label","singleOptions"});
                %colorbar(colorbar_kws_cell{:});
                %ylabel(self.currentFig.colorbar,self.colorbar_kws.label,self.colorbar_kws.fontsize);
            %else
                colorbar("off");
            end

            % add title if needed

            if ~isempty(self.title)
                title(self.title);
            end

            self.doBasePlotStuff([],[]);

        end % function plot

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef