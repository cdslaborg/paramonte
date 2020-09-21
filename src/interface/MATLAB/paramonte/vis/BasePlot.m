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
%   BasePlot(plotType,dataFrame)
%
%   This is the class for generating instances of
%   basic plots with minimally one (X)-axis. It serves as the
%   superclass for a wide variety of other multi-axes ParaMonte plots.
%
%   Parameters
%   ----------
%
%       plotType
%
%           A string indicating the name of the plot that is to be constructed.
%
%       dataFrame (optional)
%
%           A pandas dataFrame whose data will be plotted.
%
%       methodName (optional)
%
%           The name of the ParaMonte sample requesting the BasePlot.
%
%       reportEnabled (optional)
%
%           A boolean whose value indicates whether guidelines should be
%           printed in the standard output.
%
%       resetExternal (optional)
%
%           A function that resets the properties of the plot as desired
%           from outside. If provided, a pointer to this function will be
%           saved for future internal usage.
%
%   Attributes
%   ----------
%
%       dfref
%
%           A reference to the input dataFrame whose data is used to generate plots.
%           Despite its name (dfref: reference to dataFrame), this attribute is not
%           a reference to the original input dataFrame but a deep copy of it.
%           Therefore, changing its values will change the corresponding values
%           in the original dataFrame.
%
%       xcolumns
%
%           Optional property that determines the columns of dataFrame to serve as
%           the x-values. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  xcolumns = [7,8,9]
%               2.  xcolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  xcolumns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  xcolumns = 7:9      # every column in the data frame starting from column #7 to #9
%               5.  xcolumns = 7:2:20   # every other column in the data frame starting from column #7 to #20
%
%           WARNING
%
%               In all cases, xcolumns must have a length that is either 0, or 1, 
%               or equal to the length of ycolumns. If the length is 1, then xcolumns 
%               will be plotted against data corresponding to each element of ycolumns.
%               If it is an empty object having length 0, then a default value will be used.
%
%               The default value is the indices of the rows of the input dataFrame.
%
%       rows
%
%           A numeric vector that determines the rows of dataFrame
%           to be visualized. It can be either:
%
%               1.  a numeric range, or,
%               2.  a list of row indices of the dataFrame.
%
%           Example usage:
%
%               1.  rows = 15:-2:8
%               2.  rows = [12,46,7,8,9,4,7,163]
%
%           If not provided, the default includes all rows of the dataFrame.
%
%       axes
%
%           A MATLAB struct() with the following components:
%
%               kws
%
%                   A MATLAB struct() whose components will be directly passed
%                   ``set(gca,axes.kws{:})`` to set the propoerties of the
%                   current axes of the current figure. For example:
%
%           Example usage:
%
%               axes.kws.xscale = "log";
%
%           If a desired property is missing among the struct fields,
%           simply add the field and its value to axes.kws.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ``axes.kws.xscale`` and ``axes.kws.Xscale`` are the same,
%               and only one of the two will be processed.
%
%       figure
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A boolean indicating whether a call to the
%                   ``figure()`` function of MATLAB should be made or not.
%                   If a call is made, a new figure will be generated.
%                   Otherwise, the current active figure will be used.
%
%               kws
%
%                   A MATLAB struct() whose fields are directly passed to
%                   the ``figure()`` function of MATLAB.
%
%           Example usage:
%
%               figure.enabled = false; % do not make new plot
%               figure.kws.color = "none"; % set the background color to none (transparent)
%
%           If a desired property is missing among the struct fields,
%           simply add the field and its value to figure.kws.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ``axes.kws.xscale`` and ``axes.kws.Xscale`` are the same,
%               and only one of the two will be processed.
%
%       legend
%
%           A MATLAB struct() with the following components:
%
%               enabled
%
%                   A boolean indicating whether a call to the
%                   ``legend()`` function of MATLAB should be made or not.
%                   If a call is made, a new figure will be generated.
%                   Otherwise, the current active figure will be used.
%
%               kws
%
%                   A MATLAB struct() whose fields are directly passed to
%                   the ``legend()`` function of MATLAB.
%
%           Example usage:
%
%               legend.enabled = false; % do not make new plot
%               legend.kws.color = "none"; % set the background color to none (transparent)
%
%           If a desired property is missing among the struct fields,
%           simply add the field and its value to figure.kws.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ``axes.kws.xscale`` and ``axes.kws.Xscale`` are the same,
%               and only one of the two will be processed.
%
%           A MATLAB struct() whose components' values are passed to MATLAB's ``legend()`` function.
%           If your desired attribute is missing from the fieldnames of legend.kws, simply add
%           a new field named as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               legend.enabled = true; % add legend
%               legend.labels = ["this object","that object"]; % legend labels
%
%           NOTE
%
%               A legend will be added to plot only if
%               plots with no colormap are requested.
%
%           NOTE
%
%               If no legend labels is provided and legend is enabled, the names
%               of the columns of the dataFrame will be used.
%
%           WARNING
%
%               Keep in mind that MATLAB keyword arguments are case-INsensitive.
%               therefore make sure you do not add the keyword as multiple different fields.
%               For example, ``axes.kws.xscale`` and ``axes.kws.Xscale`` are the same,
%               and only one of the two will be processed.
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
%       An object of BasePlot class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef BasePlot < dynamicprops

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        dfref
        figure
        currentFig
    end

    properties (Access = protected, Hidden)
        isdryrun = [];
    end

    properties (Hidden)
        resetExternal
        rowsindex
        type
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% BasePlot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = BasePlot( plotType ...
                                , dataFrame ...
                                , resetExternal ...
                                )

            self.dfref = dataFrame;

            self.type = struct();
            self.type.name = plotType;

            plotTypeLower = lower(plotType);

            self.type.isLine       = contains(plotTypeLower,"line"      );
            self.type.isScatter    = contains(plotTypeLower,"scatter"   );
            self.type.isHeatmap    = strcmp  (plotTypeLower,"heatmap"   );
            self.type.isHistfit    = strcmp  (plotTypeLower,"histfit"   );
            self.type.isHistogram2 = strcmp  (plotTypeLower,"histogram2");
            self.type.isHistogram  = strcmp  (plotTypeLower,"histogram" );
            self.type.isEllipsoid  = contains(plotTypeLower,"covmat"    );
            self.type.isEllipsoid  = contains(plotTypeLower,"cormat"    );
            self.type.isGridPlot   = strcmp  (plotTypeLower,"grid"      );
            self.type.isContour3   = strcmp  (plotTypeLower,"contour3"  );
            self.type.isContourf   = strcmp  (plotTypeLower,"contourf"  );
            self.type.isContour    = strcmp  (plotTypeLower,"contour"   );

            self.type.is3d         = contains(plotTypeLower,"3") || self.type.isHistogram2;
            self.type.is1d         = self.type.isHistfit || self.type.isHistogram;
            self.type.is2d         = ~(self.type.isGridPlot || self.type.is1d || self.type.is3d);
            self.type.isDiffusionPlot = self.type.isContour || self.type.isContourf || self.type.isContour3;

            self.type.hasTarget = ~(self.type.isGridPlot || self.type.isHeatmap || self.type.is3d);

            %prop="rowsindex"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %if self.type.isGridPlot
            %    prop="colnames"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %    prop="colindex"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %else
            %    prop="xcolnames"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %    prop="xcolindex"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %    prop="ycolnames"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %    prop="ycolindex"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %    prop="zcolnames"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %    prop="zcolindex"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            %end

            if ~self.type.isHeatmap
                prop="rows"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="axes"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                prop="legend"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
            end

            if nargin<3 || (nargin==3 && isempty(resetExternal))
                self.resetExternal = @self.resetInternal;
            elseif nargin==3
                self.resetExternal = resetExternal;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% reset
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset  ( self ...
                        , resetType ...
                        ..., varargin ...
                        )
            %
            %   Reset the properties of the plot to the original default settings.
            %   Use this method when you change many attributes of the plot and
            %   you want to clean up and go back to the default settings.
            %
            %       Parameters
            %       ----------
            %
            %           resetType (optional)
            %
            %               An optional string with possible value of ``"hard"``.
            %               If provided, the plot object will be regenerated from scratch.
            %               This includes reading the original data frame again and resetting
            %               everything. If not provided, then only the plot settings will be
            %               reset without reseting the dataFrame.
            %
            %       Returns
            %       -------
            %
            %           None
            %
            %       Example
            %       -------
            %
            %               reset()         # reset the plot to the default settings
            %               reset("soft")   # reset the plot to the default settings
            %               reset("hard")   # regenerate the plot from scratch
            %
            if nargin<2
                resetType = "soft";
            end
            try
                self.resetExternal  ( resetType ...
                                    , self.type.name ...
                                    ..., varargin{:} ...
                                    ); % calls the external reset function
            catch
                self.resetExternal(); % calls the local reset routine
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% exportFig
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function exportFig(self, file, varargin)
            %
            %   Export the current figure to external file via the export_fig library.
            %
            %   Parameters
            %   ----------
            %
            %       file (required)
            %
            %           A string or char vector containing the path for the output generated figure file.
            %           The type of output file is determined by its extension (e.g., pdf, png, ...).
            %
            %           Example:
            %
            %               exportFig("gridplot.png")
            %
            %       varargin (optional)
            %
            %           The set of input arguments to the ``export_fig()`` function of the export_fig
            %           MATLAB library (except the path to the output exported figure which is given by `file`).
            %
            %           NOTE
            %
            %               If this input argument is missing from the input, then a default set of
            %               options, determined by the export_fig library, will be used.
            %
            %                   *.png : "-m4 -transparent"
            %                   other : []
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Example
            %   -------
            %
            %       expoortFig(); % export the figure with its name taken from the `outputFile` component of the object.
            %       expoortFig("gridplot.pdf") % export figure to a PDF file.
            %       expoortFig("gridplot.png", "-m4 -transparent") % export a large png plot of magnitude 4 with transparency.
            %
            if nargin==1 || ~( (isstring(file) || ischar(file)) && getVecLen(file) )
                error   ( "The method ``exportFig()`` takes at least one input argument ``file`` which is a " ...
                        + "string or char-vector representing the path to the output figure file, for example, " + newline ...
                        + newline ...
                        + "    expoortFig('gridplot.png')" + newline ...
                        + "    expoortFig('gridplot.png', '-m4')" + newline ...
                        + "    expoortFig('gridplot.png', '-m4 -transparent')" + newline ...
                        + newline ...
                        );
            end

            set(0, "CurrentFigure", self.currentFig.figure);
            if any(contains(string(varargin),"-transparent"))
                transparencyRequested = true;
                try
                    set(self.currentFig.figure,"color","none");
                catch
                    warning("Failed to set the color property of gcf to ""none"".");
                end
                try
                    set(gca,"color","none");
                catch
                    warning("Failed to set the color property of gca to ""none"".");
                end
            else
                transparencyRequested = false;
            end
            %for i = 1:length(varargin)
            %    if contains(varargin{i},"-transparent")
            %        transparencyRequested = true;
            %        set(self.currentFig.figure,"color","none");
            %        set(gca,"color","none");
            %    end
            %    if isa(varargin{i},"string")
            %        varargin{i} = convertStringsToChars(varargin{i});
            %    end
            %end

            export_fig(file, varargin{:});

            if transparencyRequested
                try
                    set(self.currentFig.figure,"color","default");
                catch
                    warning("failed to set the color property of gcf back to ""default"".");
                end
                try
                    set(gca,"color","default");
                catch
                    warning("failed to set the color property of gca back to ""default"".");
                end
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% getLogLinSpace
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function LogLinSpace = getLogLinSpace(self,base,logskip,lowerLim,upperLim)
            %
            %   Generate logarithmically-uniformly-spaced unique integer numbers
            %   between the input lowerLim and upperLim. These numbers are to
            %   be used as the row indices in the plots.
            %
            %   % NOTE: This function does not exist for HeatmapPlot objects.
            %
            %   Parameters
            %   ----------
            %
            %       base
            %
            %           The base of the logarithm.
            %
            %       logskip
            %
            %           The minimum logarithmic space jump between
            %           the generated log-linearly-spaced integer numbers.
            %
            %       lowerLim (optional along with upperLim)
            %
            %           The natural (non-logarithmic) lower limit of the
            %           generated log-linearly-spaced integer numbers.
            %           If not provided, the default value is 1.
            %
            %           WARNING: if lowerLim is not provided as input,
            %           WARNING: then upperLim cannot be provided either,
            %           WARNING: otherwise, the value of upperLim will be
            %           WARNING: used as the value of lowerLim and upperLim
            %           WARNING: will be assumed to not have been provided.
            %
            %       upperLim (optional)
            %
            %           The natural (non-logarithmic) upper limit of the
            %           generated log-linearly-spaced integer numbers.
            %           If not provided, the default value is the maximum
            %           of the number of the rows of the input dataframe
            %           to the BasePlot constructor.
            %
            %   Returns
            %   -------
            %
            %       A set of unique log-linearly-spaced integer numbers.
            %
            %   Example
            %   -------
            %
            %       rows = getLogLinSpace(1.01, 1, 1, 10000)
            %
            wrongSyntaxDetected = false;
            if nargin<5
                upperLim = length(self.dfref{:,1});
                if nargin<4
                    lowerLim = 1;
                    if nargin<3
                        wrongSyntaxDetected = true;
                    end
                end
            elseif nargin>5
                wrongSyntaxDetected = true;
            end
            if wrongSyntaxDetected
                error("Incorrect number of input arguments. Usage: getLogLinSpace(base,logskip,lowerLim,upperLim)")
            end
            LogLinSpace = getLogIntSpace(base,logskip,lowerLim,upperLim);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% helpme
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function helpme(self,varargin)
            if nargin==1
                object = "BasePlot";
            else
                for i = 1:nargin-1
                    if strcmpi(varargin{i},"BasePlot")
                        object = "BasePlot";
                    elseif strcmpi(varargin{i},"getLogLinSpace")
                        object = "getLogLinSpace";
                    else
                        error("The helpme() method takes only arguments with these possible values: ""BasePlot"", ""getLogLinSpace""")
                    end
                    cmd = "doc " + object;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access=public,Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetInternal(self)

            if ~(self.type.isHeatmap || self.type.isGridPlot)

                self.rows = {};

                self.legend.enabled = false;
                self.legend.labels = {};
                self.legend.kws = struct();
                self.legend.kws.box = "off";
                self.legend.kws.color = "none";
                self.legend.kws.fontSize = [];
                self.legend.kws.location = "best";
                self.legend.kws.interpreter = "none";

                self.axes.kws = struct();
                self.axes.kws.color = "white";
                self.axes.kws.xscale = "linear";
                self.axes.kws.yscale = "linear";
                self.axes.kws.zscale = "linear";

                if ~self.type.isGridPlot
                    self.axes.kws.box = "on";
                    self.axes.kws.xgrid = "on";
                    self.axes.kws.ygrid = "on";
                    if self.type.is3d
                        self.axes.kws.zgrid = "on";
                    end
                end

            end

            self.figure.enabled = true;
            self.figure.kws = struct();
            self.figure.kws.color = "white";

            self.currentFig = struct();

            if self.type.hasTarget
                prop="target"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
                self.target = Target_class();
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function parseArgs(self,varargin)

            vararginLen = length(varargin);
            selfProperties = string(properties(self));
            selfPropertiesLen = length(selfProperties);

            for i = 1:2:vararginLen % walk through input key val.
                propertyDoesNotExist = true;
                vararginItemString = string(varargin{i});
                for ip = 1:selfPropertiesLen % walk through object prop val.
                    if strcmpi(vararginItemString,selfProperties(ip))
                        propertyDoesNotExist = false;
                        if i < vararginLen % this must be here. checks for the correct pairing of key, val.
                            if isstruct(self.(selfProperties(ip))) && iscell(varargin{i+1})
                                self.(selfProperties(ip)) = parseArgs( self.(selfProperties(ip)) , varargin{i+1}{:} );
                            else
                                self.(selfProperties(ip)) = varargin{i+1};
                            end
                        else
                            error("The corresponding value for the property """ + string(selfProperties(ip)) + """ is missing as input argument.");
                        end
                        break;
                    end
                end
                if propertyDoesNotExist
                    error("The requested property """ + string(varargin{i}) + """ does not exist.");
                end
            end

        end

        % function parseArgs(self,varargin)

            % vararginLen = length(varargin);
            % for i = 1:2:vararginLen
                % propertyDoesNotExist = true;
                % selfProperties = string(properties(self));
                % selfPropertiesLen = length(selfProperties);
                % for ip = 1:selfPropertiesLen
                    % if strcmp(string(varargin{i}),string(selfProperties(ip)))
                        % propertyDoesNotExist = false;
                        % if i < vararginLen
                            % self.(selfProperties(ip)) = varargin{i+1};
                        % else
                            % error("The corresponding value for the property """ + string(selfProperties(ip)) + """ is missing as input argument.");
                        % end
                        % break;
                    % end
                % end
                % if propertyDoesNotExist
                    % error("The requested the property """ + string(varargin{i}) + """ does not exist.");
                % end
            % end

        % end % function parseArgs

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function doBasePlotStuff(self)

            if ~(self.type.isHeatmap || self.type.isGridPlot)

                % add legend

                if self.legend.enabled
                    if isa(self.legend.labels,"cell") || isa(self.legend.labels,"string")
                        %if getVecLen(self.legend.labels)~=length(lglabels)
                        %    self.legend.labels = {lglabels{:}};
                        %end
                        if isempty(self.legend.kws.fontSize)
                            self.legend.kws.fontSize = self.currentFig.xlabel.FontSize;
                        end
                        legend_kws_cell = convertStruct2Cell(self.legend.kws,{"enabled","singleOptions","labels"});
                        if self.type.isScatter && self.scatter.enabled
                            thisPlotName = "scatter"; if self.type.is3d; thisPlotName = thisPlotName + "3"; end
                            self.currentFig.legend = legend([self.currentFig.(thisPlotName){:}], self.legend.labels{:},legend_kws_cell{:});
                        elseif self.type.isHistfit && self.histfit.enabled
                            nhandle = length(self.currentFig.histfit);
                            handleCell = cell(nhandle,1);
                            for i = 1:nhandle
                                handleCell{i} = self.currentFig.histfit{i}(1); % get the handle to the histogram component of histfit
                            end
                            self.currentFig.legend = legend([handleCell{:}], self.legend.labels{:},legend_kws_cell{:});
                        else
                            self.currentFig.legend = legend(self.legend.labels{:},legend_kws_cell{:});
                        end
                        %if isfield(self.legend,"color")
                        %    set(self.currentFig.legend, 'color', self.legend.color);
                        %end
                    else
                        error   ( newline ...
                                + "The input ""legend.kws.labels"" must be a cell array of string values." ...
                                + newline ...
                                );
                    end
                else
                    legend(self.currentFig.axes,"off");
                end

                if ~isempty(self.axes.kws)
                    axes_kws_cell = convertStruct2Cell(self.axes.kws,{"enabled","singleOptions","labels"});
                    %if isfield(self.axes.kws,"singleOptions"); axes_kws_cell = { axes_kws_cell{:}, self.axes.kws.singleOptions{:} }; end
                    set(gca, axes_kws_cell{:});
                end

                if self.type.hasTarget
                    self.target.currentFig.axes = self.currentFig.axes;
                end

            end

            %%%% add title if needed

            if ~self.type.isHeatmap && (isfield(self,"title") || isprop(self,"title")) && isfield(self.title,"enabled") && self.title.enabled

                title_kws_cell = {};
                for fname = ["text", "subtext"] % do not change the order of elements here
                    if isfield(self.title,fname) && getVecLen(self.title.(fname))
                        title_kws_cell = { title_kws_cell{:}, self.title.(fname) };
                    end
                end

                if ~isempty(title_kws_cell)
                    if isfield(self.title,"kws") && ~isempty(self.title.kws) && isstruct(self.title.kws)
                        kws_cell = convertStruct2Cell(self.title.kws);
                        title_kws_cell = {title_kws_cell{:}, kws_cell{:}};
                    end
                    self.currentFig.title = title( self.currentFig.axes, title_kws_cell{:} );
                end

            end

            %if self.type.hasTarget && ~self.type.isEllipsoid
            %    xcolumnsLen = length(self.xcolumns);
            %    self.target.values = zeros(length(self.xcolumns),2);
            %    yvalues = 0;
            %    for i = 1:xcolumnsLen
            %        yvalues = 0;
            %        if ~(self.type.isHistfit || self.type.isHistogram)
            %            yvalues = mean(self.dfref{:,self.ycolindex(i)});
            %        end
            %        self.target.values(i,:) = [ mean(self.dfref{:,self.xcolindex(i)}), yvalues ];
            %    end
            %end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = private)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef BasePlot < dynamicprops
