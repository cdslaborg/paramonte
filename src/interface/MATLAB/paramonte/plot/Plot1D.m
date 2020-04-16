classdef Plot1D < dynamicprops

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)

        rows
        xcolumns
        gca_kws
        gcf_kws
        legend_kws
        currentFig
        outputFile

    end

    properties (Access = protected, Hidden)
        dfref = [];
        isdryrun = [];
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = Plot1D(varargin)

            if nargin==1
                dataFrame = varargin{1};
            else
                dataFrame = [];
            end
            self.dfref = dataFrame;
            self.reset();

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function exportFig(self,varargin)

            set(0, 'CurrentFigure', self.currentFig.gcf);
            if any(contains(string(varargin),"-transparent"))
                transparencyRequested = true;
                set(self.currentFig.gcf,'color','none');
                set(gca,'color','none');
            else
                transparencyRequested = false;
            end
            export_fig(varargin{:});
            if transparencyRequested
                set(self.currentFig.gcf,'color','default');
                set(gca,'color','default');
            end

        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

    end % methods

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function reset(self)

            self.rows = {};
            self.xcolumns = {};
            self.gca_kws = [];
            self.gcf_kws = {};
            self.legend_kws = {};
            self.currentFig = struct();
            self.outputFile = [];

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function parseArgs(self,varargin)

            vararginLen = length(varargin);
            for i = 1:2:vararginLen
                propertyDoesNotExist = true;
                selfProperties = properties(self);
                selfPropertiesLen = length(selfProperties);
                for ip = 1:selfPropertiesLen
                    if strcmp(string(varargin{i}),string(selfProperties{ip}))
                        propertyDoesNotExist = false;
                        if i < vararginLen
                            self.(selfProperties{ip}) = varargin{i+1};
                        else
                            error("The corresponding value for the property """ + string(selfProperties{ip}) + """ is missing as input argument.");
                        end
                        break;
                    end
                end
                if propertyDoesNotExist
                    error("The requested the property """ + string(varargin{i}) + """ does not exist.");
                end
            end

        end % function parseArgs

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function doStuffPlot1D(self,lgEnabled,lglabels)

            % add legend

            if lgEnabled
                if isa(self.legend_kws,"cell")
                    if isempty(self.legend_kws)
                        self.legend_kws = [ self.legend_kws , {lglabels} ];
                    end
                    [~, keyFound] = getKeyVal("fontsize",self.legend_kws{:});
                    if ~keyFound
                        fontsize_kws = {"fontsize",self.currentFig.ylabel.FontSize};
                        self.legend_kws = [ self.legend_kws , fontsize_kws ];
                    end
                    self.currentFig.legend = legend(self.legend_kws{:});
                else
                    error   ( "The input argument 'legend_kws' must be a cell array of string values." ...
                            );
                end
            else
                legend(self.currentFig.gca,'off');
            end

            if isa(self.gca_kws,"cell")
                if ~isempty(self.gca_kws)
                    set(gca, self.gca_kws{:})
                end
            end

            if isa(self.outputFile,"string") || isa(self.outputFile,"char")
                self.exportFig(self.outputFile,"-m2 -transparent");
            end

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (Access = private)

    %***************************************************************************************************************************
    %***************************************************************************************************************************

end % classdef Plot1D < dynamicprops