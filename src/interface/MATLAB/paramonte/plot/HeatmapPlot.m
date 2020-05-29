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
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HeatmapPlot(dataFrame, plotType)
%
%   This is the HeatmapPlot class for generating instances of
%   heatmap plots, built upon MATLAB's builtin function `heatmap()`.
%
%   NOTE: This is a low-level ParaMonte class and is not meant
%   NOTE: to be directly instantiated by the user.
%
%   Parameters
%   ----------
%
%       dataFrame
%
%           a MATLAB data Table from which the selected data will be plotted.
%           This is a low-level internal argument and is not meant
%           to be accessed or be provided by the user.
%
%   Attributes
%   ----------
%
%       ycolumns
%
%           optional property that determines the columns of dataFrame to serve as
%           the y-values. It can have multiple forms:
%
%               1.  a numeric or cell array of column indices in the input dataFrame.
%               2.  a string or cell array of column names in dataFrame.Properties.VariableNames.
%               3.  a cell array of a mix of the above two.
%               4.  a numeric range.
%
%           Example usage:
%
%               1.  ycolumns = [7,8,9]
%               2.  ycolumns = ["SampleLogFunc","SampleVariable1"]
%               3.  ycolumns = {"SampleLogFunc",9,"SampleVariable1"}
%               4.  ycolumns = 7:9      # every column in the data frame starting from column #7
%               5.  ycolumns = 7:2:20   # every other column in the data frame starting from column #7
%
%           WARNING: In all cases, ycolumns must have a length that is either 0, or 1, or equal
%           WARNING: to the length of xcolumns. If the length is 1, then ycolumns will be
%           WARNING: plotted against data corresponding to each element of ycolumns.
%           WARNING: If it is an empty object having length 0, then the default value will be used.
%
%           The default value is the names of all columns of the input dataFrame.
%
%       colormap
%
%           A MATLAB struct() property with two components:
%
%               1. enabled: logical value. If `true`, the colormap will be applied to the plot
%               1. values: a string or any other value that the colormap function of MATLAB accepts as input.
%
%           Example usage:
%
%               1.  colormap = "autumn"
%               1.  colormap = "winter"
%
%           If colormap is not provided or is empty, the default will be "autumn".
%
%       colorbar_kws
%
%           A MATLAB struct() whose components' values are passed to MATLAB's colorbar function.
%           If your desired attribute is missing from the fieldnames of colorbar_kws, simply add
%           a new field named as the attribute and assign the desired value to it.
%
%           Example usage:
%
%               colorbar_kws.enabled = true % add colorbar
%               colorbar_kws.location = "west"
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to colorbar_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: For example, colorbar_kws.color and colorbar_kws.Color are the same,
%           WARNING: and only one of the two will be processed.
%
%       heatmap_kws
%
%           A MATLAB struct() whose fields (with the exception of few, e.g., enabled, singleOptions, ...)
%           are directly passed to the `heatmap()` function of MATLAB.
%
%           Example usage:
%
%               heatmap_kws.enabled = true; % add heatmap()
%               heatmap_kws.ColorLimits = [-1 1];
%
%           If a desired property is missing among the struct fields, simply add the field
%           and its value to heatmap_kws.
%
%           WARNING: keep in mind that MATLAB keyword arguments are case-INsensitive.
%           WARNING: therefore make sure you do not add the keyword as multiple different fields.
%           WARNING: heatmap_kws.colorLimits and heatmap_kws.ColorLimits are the same,
%           WARNING: and only one of the two will be processed.
%
%       title
%
%           A string that is passed to the title() function of MATLAB to add title to the plot.
%           If the property is empty, not title will be added.
%
%       precision
%
%           A numeric scalar value, representing the number of digits after the decimal point 
%           for the values that appear in each cell of the heatmap.
%           The default value is not set to anything.
%
%   Superclass Attributes
%   ----------------------
%
%       See the documentation for the BasePlot class
%
%   Returns
%   -------
%
%       An object of HeatmapPlot class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
            %
            %   Adjust the limits of the colormap of the heatmap, according to 
            %   the user-input value or the default symmetric value.
            %
            %   Parameters
            %   ----------
            %
            %       newLimits (optional)
            %
            %           A vector of two scalar numeric values representing 
            %           the new limits of the colorbar. 
            %
            %           If newLimits is not given as input, the new limits will be adjusted
            %           so that the colorbar limits extend symmetrically from the  
            %           negative absolute maximum value of the current limits to the 
            %           positive absolute maximum value of the current limits.
            %
            %   Returns
            %   -------
            %
            %       None. 
            %
            %   Example
            %   -------
            %
            %       adjustColorLimits([-1 1])
            %
            if ~(isfield(self.currentFig,"heatmap") && isfield(self.currentFig.heatmap,"ColorLimits"))
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
            else
                error(newline + "There is no component ""currentFig.heatmap.ColorLimits"" for this object to adjust the limits." + newline);
            end
        end

        %***********************************************************************************************************************
        %***********************************************************************************************************************

        function plot(self,varargin)
            %
            %   Generate a plot from the selected columns of the object's dataFrame.
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
            %       plot("ycolumns",[8])
            %       plot("ycolumns",7:10)
            %

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
                %ylabel(self.currentFig.colorbar,self.colorbar_kws.label,self.colorbar_kws.fontsize, "Interpreter", "none");
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