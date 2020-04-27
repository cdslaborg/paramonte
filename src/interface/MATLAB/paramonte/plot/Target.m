classdef Target < dynamicprops
%   This is the Target class for generating instances 
%   of a target on the current active figure and axis.
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
%       value
%           a pair (list, array, or tuple) of floats, representing the (x,y) cordinates of the target.
%       axes
%           the axes object to which the target must be added.
%           The default is None in which case the output of matplotlib's gca()
%           function wil be used to get the current active axes.
%       hline_kws
%           optional dictionary of keyword arguments to be passed to matplotlib's 
%           hline() function. For example: 
%               hline = { "linewidth": 0.5,
%                           "color": "orangered",
%                           "linestyle": ":",
%                           "xmin": 0.1
%                           "xmax": 1.9
%                           }
%           The default is {}, which will be appropriately populated.
%           If set to None, no horizontal target line will be added.
%       vline_kws
%           optional dictionary of keyword arguments to be passed to matplotlib's 
%           vline() function. For example: 
%               vline = { "linewidth": 0.5,
%                           "color": "orangered",
%                           "linestyle": ":",
%                           "ymin": 0.1
%                           "ymax": 1.9
%                           }
%           The default is {}, which will be appropriately populated.
%           If set to None, no vertical target line will be added.
%       scatter_kws
%           optional dictionary of keyword arguments to be passed to matplotlib's 
%           scatter() function. For example: 
%               scatter_kws = {"s":10,"color":"green"}
%           The default is {}.
%           If set to None, no target point will be plotted.
%       outputFile
%           optional string representing the name of the output file in which 
%           the figure will be saved. If not provided, no output file will be generated.
%
%   Attributes
%   ----------
%       All of the parameters described above, except axes.
%           The input axes object (whether provided by the user or fetched by the class)
%           will be stored as a component of the object's attribute currentFig.
%       currentFig
%           an object of class ParaMonteFigure which is initially None, but upon 
%           making a plot, is populated with attributes representing all outputs 
%           from matplotlib/seaborn function calls, with the following attributes:
%               axes
%                   the output of matplotlib's gca() function or,
%                   the input axes object by the user.
%               hline
%                   the output of matplotlib's hline() function.
%               vline
%                   the output of matplotlib's vline() function.
%               scatter
%                   the output of matplotlib's scatter() function.
%
%   Returns
%   -------
%       Target
%
%   ----------------------------------------------------------------------

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        value
        hline_kws
        vline_kws
        scatter_kws
        currentFig
        xlimits
        ylimits
    end

    properties (Hidden)
        isdryrun
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function reset(self)

            xlimits = [];
            ylimits = [];
            self.value = [];
            self.hline_kws = struct();
            self.vline_kws = struct();
            self.scatter_kws = struct();

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

        function self = Target(varargin)
            self.reset()
            self.parseArgs(varargin{:});
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function plot(self,varargin)
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

            self.parseArgs(varargin{:});

            % setup the horizontal line

            %rgbColor = getRGB("deep carrot orange");
            rgbColor = getRGB("dark midnight blue");
            for fname = ["hline_kws","vline_kws","scatter_kws"]
                if contains(fname,"line")
                    if ~isfield(self.(fname),"linewidth") || isempty(self.(fname).linewidth)
                        self.(fname).linewidth = 1;
                    end
                    if ~isfield(self.(fname),"linestyle") || isempty(self.(fname).linestyle)
                        self.(fname).linestyle = "-";
                    end
                    plotName = "line()";
                else
                    if ~isfield(self.(fname),"size") || isempty(self.(fname).size)
                        self.(fname).size = 30;
                    end
                    plotName = "scatter()";
                end
                if ~isfield(self.(fname),"color") || isempty(self.(fname).color)
                    self.(fname).color = rgbColor;
                end
                if ~isfield(self.(fname),"singleOptions") || isempty(self.(fname).singleOptions)
                    self.(fname).singleOptions = {};
                elseif ~isa(self.(fname).singleOptions,"cell")
                    error   ( newline ...
                            + "The singleOptions component of " + fname + " must be a cell array of input options to the " + plotName + " function of MATLAB." ...
                            + newline ...
                            );
                end
            end
            if ~any(strcmp(self.scatter_kws.singleOptions,"filled"))
                self.scatter_kws.singleOptions = {'filled'};
            end

            % generate figure and axes if needed

            if ~isfield(self.currentFig,"gcf") || ~isgraphics(self.currentFig.gcf)
                self.currentFig.gcf = get(groot,'CurrentFigure');
            end
            self.currentFig.gca = [];
            if (~isfield(self.currentFig,"gca") || isempty(isgraphics(self.currentFig.gca))) && ~isempty(self.currentFig.gcf)
                self.currentFig.gca = self.currentFig.gcf.CurrentAxes;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "CurrentFigure", self.currentFig.gcf);
            set(self.currentFig.gcf, "CurrentAxes", self.currentFig.gca);
            hold on;

            % set what to plot

            targetExists = isa(self.value,"numeric") && length(self.value(1,:))==2;
            if ~targetExists
                error   ( newline ...
                        + "The input target value must be a pair of numeric scalars representing the X and Y coordinates of the target." ...
                        + newline ...
                        )
            end

            % check what to plot

            for fname = ["hline_kws","vline_kws","scatter_kws"]
                if ~isfield(self.(fname),"enabled")
                    self.(fname).enabled = true;
                end
                if ~isa(self.(fname).enabled,"logical")
                    error   ( newline ...
                            + "target." + fname + ".enabled must contain a logical value, true or false." ...
                            + newline ...
                            );
                end
            end

            % generate keyword options

            hline_kws_cell = {};
            fnameList = fieldnames(self.hline_kws);
            for i = 1:length(fnameList)
                if ~( strcmp(fnameList{i},"singleOptions") || strcmp(fnameList{i},"enabled") )
                    hline_kws_cell = { hline_kws_cell{:}, fnameList{i}, self.hline_kws.(fnameList{i}) };
                end
            end
            hline_kws_cell = { hline_kws_cell{:}, self.hline_kws.singleOptions{:} };

            vline_kws_cell = {};
            fnameList = fieldnames(self.vline_kws);
            for i = 1:length(fnameList)
                if ~( strcmp(fnameList{i},"singleOptions") || strcmp(fnameList{i},"enabled") )
                    vline_kws_cell = { vline_kws_cell{:}, fnameList{i}, self.vline_kws.(fnameList{i}) };
                end
            end
            vline_kws_cell = { vline_kws_cell{:}, self.vline_kws.singleOptions{:} };

            scatter_kws_cell = {};
            fnameList = fieldnames(self.scatter_kws);
            for i = 1:length(fnameList)
                if ~( strcmp(fnameList{i},"size") || strcmp(fnameList{i},"color") || strcmp(fnameList{i},"singleOptions") || strcmp(fnameList{i},"enabled") )
                    scatter_kws_cell = { scatter_kws_cell{:}, fnameList{i}, self.scatter_kws.(fnameList{i}) };
                end
            end
            scatter_kws_cell = { scatter_kws_cell{:}, self.scatter_kws.singleOptions{:} };

            % add target

            xlimCurrent = self.currentFig.gca.XLim;
            ylimCurrent = self.currentFig.gca.YLim;

            if isempty(self.xlimits); xlimits = xlimCurrent; end
            if isempty(self.ylimits); ylimits = ylimCurrent; end

            for irow = 1:length(self.value(:,1))

                if self.hline_kws.enabled
                    yline   ( ...
                            self.value(irow,2), ...
                            ... self.xlimits, [self.value(irow,2),self.value(irow,2)] ...
                            hline_kws_cell{:} ...
                            );
                end

                if self.vline_kws.enabled
                    xline   ( ...
                            self.value(irow,1), ...
                            ... [self.value(irow,1),self.value(irow,1)], self.ylimits, ...
                            vline_kws_cell{:} ...
                            );
                end

                if self.scatter_kws.enabled
                    scatter ( self.value(irow,1), self.value(irow,2) ...
                            , self.scatter_kws.size ...
                            , self.scatter_kws.color ...
                            , scatter_kws_cell{:} ...
                            );
                end

            end

            self.currentFig.gca.XLim = xlimCurrent;
            self.currentFig.gca.YLim = ylimCurrent;

            % add legend if requested

            %if legendEnabled: self.currentFig.axes.legend()

            hold off;

        end

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)

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

    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end