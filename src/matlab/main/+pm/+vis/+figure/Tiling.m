classdef Tiling < pm.vis.figure.Figure
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
    %           The output scalar object of class ``pm.vis.figure.Tiling``.
    %
    %   Interface
    %   ---------
    %
    %       plot = pm.vis.figure.Tiling(subplot);
    %
    %   Attributes
    %   ----------
    %
    %       See the list of class attributes below,
    %       also those of the superclass ``pm.vis.figure.Figure``.
    %
    properties(Access = public)
        %
        %   tiledlayout
        %
        %       A MATLAB ``struct`` whose fields and values are passed
        %       as keyword arguments to the MATLAB intrinsic ``tiledlayout()``.
        %
        tiledlayout = [];
        %
        %   subplot
        %
        %       The MATLAB cell matrix containing objects of superclass ``pm.vis.subplot.Subplot``
        %       each of which represents one subplot axes to display in the figure.
        %
        subplot = cell(0, 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %properties(Access = public, Hidden)
    %    %
    %    %       ncol
    %    %
    %    %           The MATLAB scalar whole-number whose value represents ``size(self.suplot, 2)``.
    %    %
    %    ncol = [];
    %    %
    %    %       nrow
    %    %
    %    %           The MATLAB scalar whole-number whose value represents ``size(self.suplot, 1)``.
    %    %
    %    nrow = [];
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Tiling(subplot, varargin)
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
                help("pm.vis.figure.Tiling");
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
            %       f = pm.vis.figure.Tiling.make(varargin);
            %
            %   Example
            %   -------
            %
            %       f = pm.vis.figure.Tiling();
            %       f.make()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            make@pm.vis.figure.Figure(self, varargin{:});

            %%%% Resize the figure to allow good default visualization.

            if  isempty(self.figure.position)
                maxSize = get(0, 'ScreenSize');
                figSize = self.fout.figure.Position;
                maxScale = maxSize(3:4) ./ figSize(3:4);
                newWidth = min(maxScale(2), size(self.subplot, 2)) * figSize(3) - 1;
                newHeight = min(maxScale(1), size(self.subplot, 1)) * figSize(4) - 1;
                figStart = [max(1, (maxSize(3) - newWidth) / 2), max(1, (maxSize(4) - newHeight) / 2)];
                self.fout.figure.Position = [figStart, newWidth, newHeight];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            kws = struct();
            for prop =  [ "tiledlayout" ...
                        ]
                if  isprop(self, prop)
                    kws.(prop) = self.comp2hash(prop);
                end
            end

            try
                self.fout.tiledlayout = tiledlayout(size(self.subplot, 1), size(self.subplot, 2), kws.tiledlayout{:}); % requires MATLAB R2019b.
            end
            iplt = 0;
            spinner = pm.timing.Spinner();
            for irow = 1 : size(self.subplot, 1)
                for icol = 1 : size(self.subplot, 2)
                    iplt = iplt + 1;
                    spinner.spin(iplt / numel(self.subplot)); disp([]);
                    if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.subplot.Subplot")
                        try
                            nexttile;
                        catch
                            subplot(size(self.subplot, 1), size(self.subplot, 2), iplt);
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
            %       pm.vis.figure.Tiling.reset() # reset all object properties to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self.tiledlayout.innerPosition = [];
            self.tiledlayout.outerPosition = [];
            self.tiledlayout.position = [];
            self.tiledlayout.positionConstraint = [];
            self.tiledlayout.padding = [];
            self.tiledlayout.tileSpacing = [];
            self.tiledlayout.units = [];

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
            %       This method causes side-effects by manipulating
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

            %self.nrow = size(self.subplot, 1);
            %self.ncol = size(self.subplot, 2);

            %%%% Set the default margins.

            self.setKeyVal("tiledlayout", "tileSpacing", "tight");
            %self.setKeyVal("tiledlayout", "padding", "tight");

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end