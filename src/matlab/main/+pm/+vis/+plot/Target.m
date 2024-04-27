%
%   Target(varargin)
%
%   This is the Target for generating instances of
%   targets to be added to other plots. 
%
%       NOTE
%
%           This is a low-level ParaMonte class and is not meant
%           to be directly instantiated by the user.
%
%   Parameters
%   ----------
%
%       varargin
%
%           Any attribute/value pair of Target class will be parsed and 
%           and added to the object upon construction.
%
%   Attributes
%   ----------
%
%       values
%
%           A row-vector length 2 of numeric scalars representing the X and Y coordinates of the target.
%           It can also be a matrix of shape (numTarget, 2), each row of numTarget-tows of which represents 
%           a pair of (x,y) coordinates to be added to the plot as targets.
%
%       hline
%
%           A MATLAB ``struct`` whose components' values are passed to MATLAB's `yline()` 
%           function if order to draw a horizontal line on the current active plot.
%           If your desired attribute is missing from the fieldnames of hline, simply add
%           a new field named the same as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               hline.enabled = true % add horizontal line to the current plot at the specified target value
%               hline.color = [0, 0.200000000000000, 0.400000000000000]; % line color: RGB orange
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, hline.color and hline.Color are the same, 
%           WARNING: and only one of the two will be processed.
%
%       vline
%
%           A MATLAB ``struct`` whose components' values are passed to MATLAB's `xline()` 
%           function if order to draw a vertical line on the current active plot.
%           If your desired attribute is missing from the fieldnames of vline, simply add
%           a new field named the same as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               vline.enabled = true % add vertical line to the current plot at the specified target value
%               vline.color = [0, 0.200000000000000, 0.400000000000000]; % line color: RGB orange
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, vline.color and vline.Color are the same, 
%           WARNING: and only one of the two will be processed.
%
%       scatter (available only in scatter/lineScatter/scatter3/lineScatter3 objects)
%
%           A MATLAB ``struct`` whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `scatter()` function of MATLAB.
%
%           Example usage:
%
%               scatter.enabled = true; % add scatter()
%               scatter.marker = ".";
%               scatter.size = 10;
%               scatter.color = "red";
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to scatter.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, scatter.marker and scatter.Marker are the same,
%           WARNING: and only one of the two will be processed.
%
%       xlimits
%
%           A row-vector length 2 of numeric scalars representing the limits of the hosizontal line.
%
%       ylimits
%
%           A row-vector length 2 of numeric scalars representing the limits of the vertical line.
%
%       funcout
%
%           A MATLAB ``struct`` whose fields are the outputs of various plotting tools 
%           used to make the current figure. These include the handle to the current figure (gcf), 
%           the handle to the current axes in the plot (gca), the handle to colorbar (if any exists), 
%           and other MATLAB plotting tools used to make to generate the figure.
%
%   Returns
%   -------
%
%       an object of Target class
%
%   Tip
%   ---
%
%       Question: Why did we append the name of this class with "_class"?
%       Answer: We had to do so since MATLAB has an intrinsic object named 
%       Answer: "Target" whose help and documentation overrides the help and 
%       Answer: the documentation that we have developed for this class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef Target < dynamicprops

    properties(Access = public)
        hline
        vline
        values
        scatter
        funcout
        xlimits
        ylimits
        label
    end

    properties(Hidden)
        isdryrun
        counter
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Hidden)

        function resetint(self)

            self.funcout = struct();
            self.xlimits = [];
            self.ylimits = [];
            self.values = [];

            self.counter = struct();
            self.counter.xline = 0;
            self.counter.yline = 0;
            self.counter.scatter = 0;

            self.vline = struct();
            self.hline = struct();
            self.vline.enabled = true;
            self.hline.enabled = true;

            self.scatter = struct();
            self.scatter.color = [];
            self.scatter.filled = true;
            self.scatter.enabled = true;

            self.isdryrun = true;
            self.make();
            self.isdryrun = false;

        end
        
        function hash2comp(self,varargin)
            vararginLen = length(varargin);
            for i = 1 : 2 : vararginLen
                propertyDoesNotExist = true;
                selfProperties = properties(self);
                selfPropertiesLen = length(selfProperties);
                for ip = 1 : selfPropertiesLen
                    if strcmp(string(varargin{i}), string(selfProperties{ip}))
                        propertyDoesNotExist = false;
                        if  i < vararginLen
                            self.(selfProperties{ip}) = varargin{i + 1};
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

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        function self = Target(varargin)
            self.resetint()
            self.hash2comp(varargin);
        end

        function make(self, varargin)
            %
            %   Add a target on an existing plot (the current active axes object)
            %   based on the 'values' attribute of the target object.
            %
            %   Parameters
            %   ----------
            %
            %       Any property,value pair of the object.
            %       If the property is a struct(), then its value must be given as a cell array,
            %       with consecutive elements representing the struct's property-name,property-value pairs.
            %       Note that all of these property-value pairs can be also directly set directly via the 
            %       object's attributes, before calling the make() method.
            %
            %   Returns
            %   -------
            %
            %       None. However, this method causes side-effects by manipulating
            %       the existing attributes of the object.
            %
            %   Example
            %   -------
            %
            %       make("values", [1500, 2]) % single target
            %       make("values", [1500, 2; 200, 2.5; 3500, 0.4]) % multiple targets corresponding to each row of `values`
            %       make("hline", {"linewidth", 2})
            %
            self.hash2comp(varargin);
            rgbColor = pm.vis.color.rgb("deep carrot orange");
            %rgbColor = pm.vis.color.rgb("dark midnight blue");

            %%%% setup the horizontal line

            for fname = ["hline", "vline", "scatter"]

                if contains(fname, "line")

                    plotName = "line()";

                    key = "linewidth"; val = 1.5; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "lineStyle"; val = "-"; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "color"; val = rgbColor; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end

                elseif contains(fname, "scatter")

                    plotName = "scatter()";

                    key = "marker"; val = "o"; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "filled"; val = true; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "color"; val = rgbColor; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "size"; val = 30; if ~isfield(self.(fname), key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end

                end

                %if ~isfield(self.(fname), "singleOptions") || isempty(self.(fname).singleOptions)
                %    self.(fname).singleOptions = {};
                %elseif ~isa(self.(fname).singleOptions, "cell")
                %    error   ( newline ...
                %            + "The singleOptions component of " + fname + " must be a cell array of input options to the " + plotName + " function of MATLAB." ...
                %            + newline ...
                %            );
                %end

            end

            %if ~any(strcmp(self.scatter.singleOptions, "filled"))
            %    self.scatter.singleOptions = {'filled'};
            %end

            %%%% generate figure and axes if needed

            if ~isfield(self.funcout, "gcf") || ~all(isgraphics(self.funcout.figure)) || isempty(isgraphics(self.funcout.figure))
                self.funcout.figure = get(groot,'CurrentFigure');
            end

            if (~isfield(self.funcout, "axes") || isempty(isgraphics(self.funcout.axes))) && ~isempty(self.funcout.figure)
                self.funcout.axes = self.funcout.figure.CurrentAxes;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            noFigureDetected = false;
            try 
                set(0, "CurrentFigure", self.funcout.figure);
            catch
                noFigureDetected = true;
            end

            %%%% if there is really no figure, report error

            if noFigureDetected || ~all(isgraphics(self.funcout.figure))
                error(newline + "There is no figure to which the target could be added. Please make a plot first." + newline);
            end

            set(self.funcout.figure, "CurrentAxes", self.funcout.axes);
            hold on;

            %%%% set what to plot

            targetExists = isa(self.values, "numeric") && length(self.values(1,:))==2;
            if ~targetExists
                error   ( newline ...
                        + "The input target values must be a pair of numeric scalars representing the X and Y coordinates of the target, " ...
                        + "or a matrix of shape (numTarget, 2), each row of numTarget-tows of which is a pair of (x,y) coordinates " ...
                        + "to be added as targets to the plot." ...
                        + newline ...
                        )
            end

            %%%% check what to plot

            for fname = ["hline", "vline", "scatter"]
                if ~isfield(self.(fname), "enabled")
                    self.(fname).enabled = true;
                end
                if ~isa(self.(fname).enabled, "logical")
                    error   ( newline ...
                            + "target." + fname + ".enabled must contain a logical value, true or false." ...
                            + newline ...
                            );
                end
            end

            %%%% generate hline keyword options

            hline_kws_cell = {};
            fnameList = fieldnames(self.hline);
            for i = 1 : length(fnameList)
                if ~strcmpi(fnameList{i}, "enabled")
                    hline_kws_cell = {hline_kws_cell{:}, fnameList{i}, self.hline.(fnameList{i})};
                end
            end

            %%%% generate vline keyword options

            vline_kws_cell = {};
            fnameList = fieldnames(self.vline);
            for i = 1 : length(fnameList)
                if ~strcmpi(fnameList{i}, "enabled")
                    vline_kws_cell = {vline_kws_cell{:}, fnameList{i}, self.vline.(fnameList{i})};
                end
            end

            %%%% generate scatter keyword options

            scatter_kws_cell = {};
            if  self.scatter.filled
                scatter_kws_cell = {"filled"};
            end
            fnameList = fieldnames(self.scatter);
            for i = 1 : length(fnameList)
                %if ~(strcmp(fnameList{i}, "size") || strcmp(fnameList{i}, "color") || strcmp(fnameList{i}, "singleOptions") || strcmp(fnameList{i}, "enabled"))
                scatter_kws_cell = {scatter_kws_cell{:}, fnameList{i}, self.scatter.(fnameList{i})};
                %end
            end

            %%%% add target

            xlimCurrent = self.funcout.axes.XLim;
            ylimCurrent = self.funcout.axes.YLim;

            if isempty(self.xlimits); xlimits = xlimCurrent; end
            if isempty(self.ylimits); ylimits = ylimCurrent; end

            for irow = 1 : length(self.values(:, 1))

                isnanx = isnan(self.values(irow, 1));
                isnany = isnan(self.values(irow, 2));

                if  self.hline.enabled && ~isnany
                    self.counter.yline = self.counter.yline + 1;
                    self.funcout.yline{self.counter.yline} = yline   ( self.values(irow, 2) ...
                                                                        ..., self.xlimits, [self.values(irow, 2), self.values(irow, 2)] ...
                                                                        , hline_kws_cell{:} ...
                                                                        );
                end

                if  self.vline.enabled && ~isnanx
                    self.counter.xline = self.counter.xline + 1;
                    self.funcout.xline{self.counter.xline} = xline   ( self.values(irow, 1) ...
                                                                        ..., [self.values(irow, 1), self.values(irow, 1)], self.ylimits ...
                                                                        , vline_kws_cell{:} ...
                                                                        );
                end

                if  self.scatter.enabled && ~(isnanx || isnany)
                    self.counter.scatter = self.counter.scatter + 1;
                    self.funcout.scatter{self.counter.scatter} = scatter ( self.values(irow, 1), self.values(irow, 2) ...
                                                                            , self.scatter.size ...
                                                                            , self.scatter.color ...
                                                                            , self.scatter.marker ...
                                                                            , scatter_kws_cell{:} ...
                                                                            );
                end

            end

            self.funcout.axes.XLim = xlimCurrent;
            self.funcout.axes.YLim = ylimCurrent;

            % add legend if requested

            %if legendEnabled: self.funcout.axes.legend()

            hold off;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end