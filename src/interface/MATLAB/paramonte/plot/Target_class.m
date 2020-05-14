%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see,
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms,
%   if you use any parts of this library for any purposes,
%   we ask you to acknowledge the ParaMonte library's usage
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Target_class(varargin)
%
%   This is the Target_class for generating instances of
%   targets to be added to other plots. 
%
%   NOTE: This is a low-level ParaMonte class and is not meant
%   NOTE: to be directly instantiated by the user.
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
%       value
%
%           A row-vector length 2 of numeric scalars representing the X and Y coordinates of the target.
%           It can also be a matrix of shape (numTarget,2), each row of numTarget-tows of which represents 
%           a pair of (x,y) coordinates to be added to the plot as targets.
%
%       hline_kws
%
%           A MATLAB struct() whose components' values are passed to MATLAB's `yline()` 
%           function if order to draw a horizontal line on the current active plot.
%           If your desired attribute is missing from the fieldnames of hline_kws, simply add
%           a new field named the same as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               hline_kws.enabled = true % add horizontal line to the current plot at the specified target value
%               hline_kws.color = [0, 0.200000000000000, 0.400000000000000]; % line color: RGB orange
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, hline_kws.color and hline_kws.Color are the same, 
%           WARNING: and only one of the two will be processed.
%
%       vline_kws
%
%           A MATLAB struct() whose components' values are passed to MATLAB's `xline()` 
%           function if order to draw a vertical line on the current active plot.
%           If your desired attribute is missing from the fieldnames of vline_kws, simply add
%           a new field named the same as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               vline_kws.enabled = true % add vertical line to the current plot at the specified target value
%               vline_kws.color = [0, 0.200000000000000, 0.400000000000000]; % line color: RGB orange
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, vline_kws.color and vline_kws.Color are the same, 
%           WARNING: and only one of the two will be processed.
%
%       scatter_kws (available only in scatter/lineScatter/scatter3/lineScatter3 objects)
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `scatter()` function of MATLAB.
%
%           Example usage:
%
%               scatter_kws.enabled = true; % add scatter()
%               scatter_kws.marker = ".";
%               scatter_kws.size = 10;
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to scatter_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, scatter_kws.marker and scatter_kws.Marker are the same,
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
%       Answer: documentation that we have developed for this class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef Target_class < dynamicprops

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Target_class(varargin)
            self.reset()
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

        function plot(self,varargin)
            %
            %   Add a target on an existing plot (the current active axes object)
            %   based on the 'value' attribute of the target object.
            %
            %   Parameters
            %   ----------
            %
            %       Any property,value pair of the object.
            %       If the property is a struct(), then its value must be given as a cell array,
            %       with consecutive elements representing the struct's property-name,property-value pairs.
            %       Note that all of these property-value pairs can be also directly set directly via the 
            %       object's attributes, before calling the plot() method.
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
            %       plot("value",[1500,2]) % single target
            %       plot("value",[1500,2;200,2.5;3500,0.4]) % multiple targets corresponding to each row of `value`
            %       plot("hline_kws",{"linewidth",2})
            %
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

            if ~isfield(self.currentFig,"gcf") || ~all(isgraphics(self.currentFig.gcf)) || isempty(isgraphics(self.currentFig.gcf))
                self.currentFig.gcf = get(groot,'CurrentFigure');
            end

            self.currentFig.gca = [];
            if (~isfield(self.currentFig,"gca") || isempty(isgraphics(self.currentFig.gca))) && ~isempty(self.currentFig.gcf)
                self.currentFig.gca = self.currentFig.gcf.CurrentAxes;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.isdryrun; return; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            noFigureDetected = false;
            try 
                set(0, "CurrentFigure", self.currentFig.gcf);
            catch
                noFigureDetected = true;
            end

            % if there is really no figure, report error

            if noFigureDetected || ~all(isgraphics(self.currentFig.gcf))
                error(newline + "There is no figure to which the target could be added. Please make a plot first." + newline);
            end

            set(self.currentFig.gcf, "CurrentAxes", self.currentFig.gca);
            hold on;

            % set what to plot

            targetExists = isa(self.value,"numeric") && length(self.value(1,:))==2;
            if ~targetExists
                error   ( newline ...
                        + "The input target value must be a pair of numeric scalars representing the X and Y coordinates of the target, " ...
                        + "or a matrix of shape (numTarget,2), each row of numTarget-tows of which is a pair of (x,y) coordinates " ...
                        + "to be added as targets to the plot." ...
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