%>  \brief
%>  This is the base class for generating instances
%>  of figures containing a tile of subplots.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.vis.figure.Figure](@ref Figure).<br>
%>
%>  \note
%>  See also the documentation of the constructor of the class [pm.vis.Tiling::Tiling](@ref Tiling::Tiling).<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 9:20 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
classdef Tiling < pm.vis.figure.Figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``tiledlayout``
        %>
        %>  A MATLAB ``struct`` whose fields and values are passed
        %>  as keyword arguments to the MATLAB intrinsic ``tiledlayout()``.<br>
        %>
        tiledlayout = [];
        %>
        %>  ``subplot``
        %>
        %>  The MATLAB cell matrix containing objects of superclass [pm.vis.Subplot](@ref Subplot)
        %>  each of which represents one subplot axes to display in the figure.<br>
        %>
        subplot = cell(0, 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Construct and return an object of class [pm.vis.Tiling](@ref Tiling).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.vis.Tiling](@ref Tiling).<br>
        %>
        %>  \param[in]  subplot     :   The input cell matrix of MATLAB objects of superclass [pm.vis.Subplot](@ref Subplot).<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.Tiling](@ref Tiling).<br>
        %>
        %>  \interface{Tiling}
        %>  \code{.m}
        %>
        %>      t = pm.vis.Tiling(subplot);
        %>      t = pm.vis.Tiling(subplot, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the list of class attributes below,
        %>  also those of the superclass [pm.vis.figure.Figure](@ref Figure).<br>
        %>
        %>  \final{Tiling}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:20 AM, University of Texas at Arlington<br>
        function self = Tiling(subplot, varargin)
            if  nargin < 1
                subplot = cell(0, 0);
            end
            failed = ~iscell(subplot);
            if ~failed
                for irow = 1 : size(subplot, 1)
                    for icol = 1 : size(subplot, 2)
                        failed = ~isempty(subplot{irow, icol}) && ~pm.introspection.istype(subplot{irow, icol}, "pm.vis.Subplot");
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
                help("pm.vis.Tiling");
                error   ( newline ...
                        + "The input argument ``subplot`` must be a MATLAB cell matrix of " + newline ...
                        + "empty objects or objects of superclass [pm.vis.Subplot](@ref Subplot)." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            self = self@pm.vis.figure.Figure(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the figure settings and specifications,
        %>  make the figure and the subplots, and return nothing.<br>
        %>
        %>  \details
        %>  The subplots are made by calling their ``make()`` methods.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      t = pm.vis.Tiling(subplot, varargin);
        %>      t.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:24 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function make(self, varargin)

            make@pm.vis.figure.Figure(self, varargin{:});

            %%%% Resize the figure to allow good default visualization.

            if  isempty(self.figure.outerPosition) && 0 < numel(self.subplot)
                fullSize = get(0, 'ScreenSize');
                maxSize = fullSize;
                maxSize(1:2) = .05 * maxSize(3:4);
                maxSize(3:4) = maxSize(3:4) - maxSize(1:2);
                figSize = self.fout.figure.OuterPosition;
                maxScale = maxSize(3:4) ./ figSize(3:4);
                newWidth = figSize(3) * min(maxScale(1), size(self.subplot, 2));
                newHeight = figSize(4) * min(maxScale(2), size(self.subplot, 1));
                figStart = [(fullSize(3) - newWidth) / 2, (fullSize(4) - newHeight) / 2];
                set(self.fout.figure, "OuterPosition", [figStart(1:2), newWidth, newHeight]);
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
            timer = pm.timing.Timer();
            spinner = pm.timing.Spinner();
            for irow = 1 : size(self.subplot, 1)
                for icol = 1 : size(self.subplot, 2)
                    iplt = iplt + 1;
                    spinner.spin(iplt / numel(self.subplot));
                    if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.Subplot")
                        try
                            nexttile;
                        catch
                            subplot(size(self.subplot, 1), size(self.subplot, 2), iplt);
                        end
                        self.subplot{irow, icol}.make();
                    end
                end
            end
            if  0 < iplt
                disp("done in " + sprintf("%.6f", string(timer.toc())) + " seconds.");
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the figure to the original default settings.<br>
        %>
        %>  \details
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.<br>
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      t = pm.vis.Tiling(subplot, varargin)
        %>      t.reset(varargin); % reset all object properties to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:25 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function reset(self, varargin)

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

        %>  \brief
        %>  Preset the tiling settings before making it.<br>
        %>
        %>  \warning
        %>  This method causes side-effects by manipulating
        %>  the existing attributes of the object.<br>
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.<br>
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      t = pm.vis.Tiling(subplot, varargin);
        %>      t.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 9:28 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, July 7 2024, 12:53 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function premake(self, varargin)

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