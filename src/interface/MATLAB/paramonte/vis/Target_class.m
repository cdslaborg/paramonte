%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Target_class(varargin)
%
%   This is the Target_class for generating instances of
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
%           It can also be a matrix of shape (numTarget,2), each row of numTarget-tows of which represents 
%           a pair of (x,y) coordinates to be added to the plot as targets.
%
%       hline
%
%           A MATLAB struct() whose components' values are passed to MATLAB's `yline()` 
%           function if order to draw a horizontal line on the current active plot.
%           If your desired attribute is missing from the fieldnames of hline.kws, simply add
%           a new field named the same as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               hline.enabled = true % add horizontal line to the current plot at the specified target value
%               hline.kws.color = [0, 0.200000000000000, 0.400000000000000]; % line color: RGB orange
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, hline.kws.color and hline.kws.Color are the same, 
%           WARNING: and only one of the two will be processed.
%
%       vline.kws
%
%           A MATLAB struct() whose components' values are passed to MATLAB's `xline()` 
%           function if order to draw a vertical line on the current active plot.
%           If your desired attribute is missing from the fieldnames of vline.kws, simply add
%           a new field named the same as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               vline.enabled = true % add vertical line to the current plot at the specified target value
%               vline.kws.color = [0, 0.200000000000000, 0.400000000000000]; % line color: RGB orange
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, vline.kws.color and vline.kws.Color are the same, 
%           WARNING: and only one of the two will be processed.
%
%       scatter.kws (available only in scatter/lineScatter/scatter3/lineScatter3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `scatter()` function of MATLAB.
%
%           Example usage:
%
%               scatter.enabled = true; % add scatter()
%               scatter.kws.marker = ".";
%               scatter.size = 10;
%               scatter.color = "red";
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to scatter.kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, scatter.kws.marker and scatter.kws.Marker are the same,
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
%       currentFig
%
%           A MATLAB struct() whose fields are the outputs of various plotting tools 
%           used to make the current figure. These include the handle to the current figure (gcf), 
%           the handle to the current axes in the plot (gca), the handle to colorbar (if any exists), 
%           and other MATLAB plotting tools used to make to generate the figure.
%
%   Returns
%   -------
%
%       an object of Target_class class
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
classdef Target_class < dynamicprops

    properties (Access = public)
        hline
        vline
        values
        scatter
        currentFig
        xlimits
        ylimits
        label
    end

    properties (Hidden)
        isdryrun
        counter
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetInternal(self)

            self.xlimits = [];
            self.ylimits = [];
            self.values = [];

            self.vline = struct();
            self.hline = struct();
            self.scatter = struct();

            self.counter = struct();
            self.counter.xline = 0;
            self.counter.yline = 0;
            self.counter.scatter = 0;
            self.currentFig = struct();

            self.vline.enabled = true;
            self.hline.enabled = true;
            self.scatter.enabled = true;

            self.vline.kws = struct();
            self.hline.kws = struct();
            self.scatter.kws = struct();
            self.scatter.color = [];
            self.scatter.filled = true;

            self.isdryrun = true;
            self.make();
            self.isdryrun = false;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Target_class(varargin)
            self.resetInternal()
            self.parseArgs(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function helpme(self,varargin)
            %
            %   Open the documentation for the input object's name in string format, otherwise, 
            %   open the documentation page for the class of the object owning the helpme() method.
            %
            %   Parameters
            %   ----------
            %
            %       This function takes at most one string argument, 
            %       which is the name of the object for which help is needed.
            %
            %   Returns
            %   -------
            %
            %       None. 
            %
            %   Example
            %   -------
            %
            %       helpme("plot")
            %
            methodNotFound = true;
            if nargin==2
                if strcmpi(varargin{1},"target")
                    cmd = "doc Target_class";
                    methodNotFound = false;
                else
                    methodList = ["plot","helpme"];
                    for method = methodList
                        if strcmpi(varargin{1},method)
                            methodNotFound = false;
                            cmd = "doc self." + method;
                        end
                    end
                end
            elseif nargin~=1
                error("The helpme() method takes at most one argument that must be string.");
            end
            if methodNotFound
                cmd = "doc self";
            end
            eval(cmd);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self,varargin)
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
            %       make("values",[1500,2]) % single target
            %       make("values",[1500,2;200,2.5;3500,0.4]) % multiple targets corresponding to each row of `values`
            %       make("hline.kws",{"lineWidth",2})
            %

            self.parseArgs(varargin{:});
            rgbColor = getRGB("deep carrot orange");
            %rgbColor = getRGB("dark midnight blue");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% setup the horizontal line
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for fname = ["hline","vline","scatter"]

                if contains(fname,"line")

                    plotName = "line()";

                    key = "lineWidth"; val = 1.5; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                    key = "lineStyle"; val = "-"; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end
                    key = "color"; val = rgbColor; if ~isfield(self.(fname).kws,key) || isempty(self.(fname).kws.(key)); self.(fname).kws.(key) = val; end

                elseif contains(fname,"scatter")

                    plotName = "scatter()";

                    key = "marker"; val = "o"; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "filled"; val = true; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "color"; val = rgbColor; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end
                    key = "size"; val = 30; if ~isfield(self.(fname),key) || isempty(self.(fname).(key)); self.(fname).(key) = val; end

                end

                %if ~isfield(self.(fname),"singleOptions") || isempty(self.(fname).singleOptions)
                %    self.(fname).singleOptions = {};
                %elseif ~isa(self.(fname).singleOptions,"cell")
                %    error   ( newline ...
                %            + "The singleOptions component of " + fname + " must be a cell array of input options to the " + plotName + " function of MATLAB." ...
                %            + newline ...
                %            );
                %end

            end

            %if ~any(strcmp(self.scatter.singleOptions,"filled"))
            %    self.scatter.singleOptions = {'filled'};
            %end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate figure and axes if needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ~isfield(self.currentFig,"gcf") || ~all(isgraphics(self.currentFig.figure)) || isempty(isgraphics(self.currentFig.figure))
                self.currentFig.figure = get(groot,'CurrentFigure');
            end

            if (~isfield(self.currentFig,"axes") || isempty(isgraphics(self.currentFig.axes))) && ~isempty(self.currentFig.figure)
                self.currentFig.axes = self.currentFig.figure.CurrentAxes;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            noFigureDetected = false;
            try 
                set(0, "CurrentFigure", self.currentFig.figure);
            catch
                noFigureDetected = true;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% if there is really no figure, report error
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if noFigureDetected || ~all(isgraphics(self.currentFig.figure))
                error(newline + "There is no figure to which the target could be added. Please make a plot first." + newline);
            end

            set(self.currentFig.figure, "CurrentAxes", self.currentFig.axes);
            hold on;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            targetExists = isa(self.values,"numeric") && length(self.values(1,:))==2;
            if ~targetExists
                error   ( newline ...
                        + "The input target values must be a pair of numeric scalars representing the X and Y coordinates of the target, " ...
                        + "or a matrix of shape (numTarget,2), each row of numTarget-tows of which is a pair of (x,y) coordinates " ...
                        + "to be added as targets to the plot." ...
                        + newline ...
                        )
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% check what to plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for fname = ["hline","vline","scatter"]
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate hline keyword options
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            hline_kws_cell = {};
            fnameList = fieldnames(self.hline.kws);
            for i = 1:length(fnameList)
                if ~( strcmp(fnameList{i},"singleOptions") || strcmp(fnameList{i},"enabled") )
                    hline_kws_cell = { hline_kws_cell{:}, fnameList{i}, self.hline.kws.(fnameList{i}) };
                end
            end
            %hline_kws_cell = { hline_kws_cell{:}, self.hline.singleOptions{:} };

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate vline keyword options
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            vline_kws_cell = {};
            fnameList = fieldnames(self.vline.kws);
            for i = 1:length(fnameList)
                if ~( strcmp(fnameList{i},"singleOptions") || strcmp(fnameList{i},"enabled") )
                    vline_kws_cell = { vline_kws_cell{:}, fnameList{i}, self.vline.kws.(fnameList{i}) };
                end
            end
            %vline_kws_cell = { vline_kws_cell{:}, self.vline.singleOptions{:} };

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% generate scatter keyword options
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            scatter_kws_cell = {};
            if self.scatter.filled; scatter_kws_cell = {"filled"}; end
            fnameList = fieldnames(self.scatter.kws);
            for i = 1:length(fnameList)
                %if ~( strcmp(fnameList{i},"size") || strcmp(fnameList{i},"color") || strcmp(fnameList{i},"singleOptions") || strcmp(fnameList{i},"enabled") )
                    scatter_kws_cell = { scatter_kws_cell{:}, fnameList{i}, self.scatter.kws.(fnameList{i}) };
                %end
            end
            %scatter_kws_cell = { scatter_kws_cell{:}, self.scatter.singleOptions{:} };

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% add target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            xlimCurrent = self.currentFig.axes.XLim;
            ylimCurrent = self.currentFig.axes.YLim;

            if isempty(self.xlimits); xlimits = xlimCurrent; end
            if isempty(self.ylimits); ylimits = ylimCurrent; end

            for irow = 1:length(self.values(:,1))

                isnanx = isnan(self.values(irow,1));
                isnany = isnan(self.values(irow,2));

                if self.hline.enabled && ~isnany
                    self.counter.yline = self.counter.yline + 1;
                    self.currentFig.yline{self.counter.yline} = yline   ( ...
                                                                        self.values(irow,2), ...
                                                                        ... self.xlimits, [self.values(irow,2),self.values(irow,2)] ...
                                                                        hline_kws_cell{:} ...
                                                                        );
                end

                if self.vline.enabled && ~isnanx
                    self.counter.xline = self.counter.xline + 1;
                    self.currentFig.xline{self.counter.xline} = xline   ( ...
                                                                        self.values(irow,1), ...
                                                                        ... [self.values(irow,1),self.values(irow,1)], self.ylimits, ...
                                                                        vline_kws_cell{:} ...
                                                                        );
                end

                if self.scatter.enabled && ~(isnanx || isnany)
                    self.counter.scatter = self.counter.scatter + 1;
                    self.currentFig.scatter{self.counter.scatter} = scatter ( self.values(irow,1), self.values(irow,2) ...
                                                                            , self.scatter.size ...
                                                                            , self.scatter.color ...
                                                                            , self.scatter.marker ...
                                                                            , scatter_kws_cell{:} ...
                                                                            );
                end

            end

            self.currentFig.axes.XLim = xlimCurrent;
            self.currentFig.axes.YLim = ylimCurrent;

            % add legend if requested

            %if legendEnabled: self.currentFig.axes.legend()

            hold off;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end