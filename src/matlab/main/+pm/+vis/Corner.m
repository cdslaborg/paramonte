classdef Corner < pm.vis.figure.Figure
    %
    %   This is the abstract class for generating instances
    %   of figures containing a tile of subplots.
    %
    %   Parameters
    %   ----------
    %
    %       subplot
    %
    %           The input cell matrix of MATLAB objects of superclass ``pm.vis.subplot.Subplot``.
    %
    %       varargin
    %
    %           Any ``property, value`` pair of the parent object.
    %           If the property is a ``struct()``, then its value must be given as a cell array,
    %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
    %           Note that all of these property-value pairs can be also directly set via the
    %           parent object attributes, before calling the ``make()`` method.
    %
    %   Returns
    %   -------
    %
    %       self
    %
    %           The output scalar object of class ``pm.vis.Corner``.
    %
    %   Interface
    %   ---------
    %
    %       plot = pm.vis.Corner(subplot);
    %
    %   Attributes
    %   ----------
    %
    %       See the list of class attributes below,
    %       also those of the superclass ``pm.vis.figure.Figure``.
    %
    properties(Access = public)
        %
        %       margin
        %
        %           The MATLAB ``struct`` whose components contain information about the figure subplot margins.
        %           This information is used to place the input subplots to the constructor of the object
        %           in the appropriate locations in the figure.
        %
        %           \warning
        %
        %               All specified margin values must be a number between zero and one.
        %               Furthermore, the specified values must make sense and be reasonable.
        %
        margin = struct("spacex", [], "spacey", [], "bottom", [], "right", [], "left", [], "top", []);
        %
        %       subplot
        %
        %           The MATLAB cell matrix containing objects of superclass ``pm.vis.subplot.Subplot``
        %           each of which represents one subplot axes to display in the figure.
        %
        subplot = cell(0, 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public, Hidden)
        %
        %       ncol
        %
        %           The MATLAB scalar whole-number whose value represents ``size(self.suplot, 2)``.
        %
        ncol = [];
        %
        %       nrow
        %
        %           The MATLAB scalar whole-number whose value represents ``size(self.suplot, 1)``.
        %
        nrow = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Corner(subplot, varargin)
            if  nargin < 1
                subplot = cell(0, 0);
            end
            failed = ~iscell(subplot);
            if ~failed
                for irow = 1 : size(subplot, 1)
                    for icol = 1 : size(subplot, 2)
                        failed = ~isempty(subplot{irow, icol}) && ~pm.introspection.istype(subplot{irow, icol}, "pm.vis.subplot.Subplot");
                        if  failed
                            break;
                        end
                    end
                    if  failed
                        break;
                    end
                end
            end
            if ~failed
                varargin = {"subplot", subplot, varargin{:}};
            else
                help("pm.vis.Corner");
                error   ( newline ...
                        + "The input argument ``subplot`` must be a MATLAB cell matrix of " + newline ...
                        + "empty objects or objects of superclass ``pm.vis.subplot.Subplot``." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            self = self@pm.vis.figure.Figure(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
            %
            %   Configure the figure settings and specifications,
            %   make the figure and the subplots, and return nothing.
            %   The subplots are made by calling their ``make()`` methods.
            %
            %   \warning
            %
            %       This method has side-effects by manipulating
            %       the existing attributes of the parent object.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes, before calling the ``make()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       f = pm.vis.Corner.make(varargin);
            %
            %   Example
            %   -------
            %
            %       f = pm.vis.Corner();
            %       f.make()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            make@pm.vis.figure.Figure(self, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            iplot = 0;
            try
                notiles = false;
                tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'none'); % requires MATLAB R2019b.
            end
            for irow = 1 : self.nrow
                for icol = 1 : self.ncol
                    iplot = iplot + 1;
                    if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.subplot.Subplot")
                        try
                            nexttile
                        catch
                            subplot(self.nrow, self.ncol, iplot);
                        end
                        self.subplot{irow, icol}.make();
                    end
                end
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self, varargin)
            %
            %   Reset the properties of the figure to the original default settings.
            %   Use this method when you change many attributes of the plot and
            %   you want to clean up and go back to the default settings.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes, before calling the ``make()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.Corner.reset() # reset all object properties to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            reset@pm.vis.figure.Figure(self, varargin{:}); % this will automatically call the ``premake()`` method of the object.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the make() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.premake(varargin{:}); % This is the subclass method!
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function premake(self, varargin)
            %
            %   Update the layout of the tile with the new axes position changes.
            %
            %   \warning
            %
            %       this method causes side-effects by manipulating
            %       the existing attributes of the object.
            %
            %   Parameters
            %   ----------
            %
            %       varargin
            %
            %           Any ``property, value`` pair of the parent object.
            %           If the property is a ``struct()``, then its value must be given as a cell array,
            %           with consecutive elements representing the struct ``property-name, property-value`` pairs.
            %           Note that all of these property-value pairs can be also directly set via the
            %           parent object attributes, before calling the ``make()`` method.
            %
            %   Returns
            %   -------
            %
            %       None.
            %
            %   Example
            %   -------
            %
            %       For example, change the left/bottom margin of the main
            %       axis of the figure to provide room for lengthy variable names.
            %       Then call the ``self.update()`` method to reflect the changes.
            %
            premake@pm.vis.figure.Figure(self, varargin{:});

            self.nrow = size(self.subplot, 1);
            self.ncol = size(self.subplot, 2);

            %%%% Set the default margins.

            if  isempty(self.margin.spacex)
                self.margin.spacex = 0.015; % x-space between subplots
            end
            if  isempty(self.margin.spacey)
                self.margin.spacey = 0.015; % y-space between subplots
            end
            if  isempty(self.margin.bottom)
                self.margin.bottom = 0.07; % 0.14;
            end
            if  isempty(self.margin.right)
                self.margin.right = 0.10;
            end
            if  isempty(self.margin.left)
                self.margin.left = 0.07; %0.1;
            end
            if  isempty(self.margin.top)
                self.margin.top = 0.0;
            end

            %for irow = 1 : self.nrow
            %    for icol = 1 : self.ncol
            %        jrow = self.nrow - irow + 1;
            %        if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.subplot.Subplot")
            %            if ~isfield(self.subplot{irow, icol}.axes, "position") || isempty(self.subplot{irow, icol}.axes.position)
            %                self.subplot{irow, icol}.axes.position = self.getpos(jrow, icol);
            %            end
            %        end
            %    end
            %end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function position = getpos(self, irow, icol)
            if  nargin == 3
                width  = max(0, (1 - self.margin.left - self.margin.right) / self.nrow - self.margin.spacex);
                height = max(0, (1 - self.margin.top - self.margin.bottom) / self.ncol - self.margin.spacey);
                position =  [ (icol - 1) * (self.margin.spacex + width)  + self.margin.left ...
                            , (irow - 1) * (self.margin.spacey + height) + self.margin.bottom ...
                            , width ...
                            , height ...
                            ];
            else
                position =  [ self.margin.right ...
                            , self.margin.bottom ...
                            , 1 - self.margin.right ...width
                            , 1 - self.margin.top ... height
                            ];
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end