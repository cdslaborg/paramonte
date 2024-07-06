%>  \brief
%>  This is the abstract class for generating instances of objects
%>  that contain the specifications of various types of tile figures.<br>
%>  This is a generic class for generating figures
%>  containing multiple subplots (axes) of the same class.
classdef Tile < [pm.vis.figure.Tiling](@ref Tiling)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  \param  template    :   The scalar object of superclass [pm.vis.subplot.Subplot](@ref Subplot)
        %>                          representing the template of the set of subplots to display.
        %>
        template = [];
        %>
        %>  \param  tileshape   :   The MATLAB vector ``[nrow, ncol]`` representing
        %>                          the number of rows and columns in the tiled layout of subplots.
        %>                          The default is the closest values of ``nrow`` and ``ncol`` whose
        %>                          multiplication yields the maximum number of data columns in the
        %>                          visualization.
        %>
        tileshape = [];
        %>
        %>  \param  tileindex   :   The MATLAB vector of size ``0`` or ``prod(tileshape)`` containing
        %>                          the list of indices of the tiling to fill with subplots, starting
        %>
        %tileindex = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>
        %>  \note
        %>  This class is not meant to be used directly by the end users.
        %>  Instead, use the subclasses of this abstract class for visualizations.
        %>
        %>  \param[in]  template    :   The input scalar object of superclass [pm.vis.subplot.Subplot](@ref Subplot).
        %>                              It serves as the template based upon which all subplots are constructed.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``template`` component of the parent object.
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class [pm.vis.tile.Tile](@ref Tile).
        %>
        %>  \interface{Tile}
        %>  \code{.m}
        %>
        %>      tile = pm.vis.tile.Tile(template);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the list of class attributes below,
        %>  also those of the superclass [pm.vis.figure.Tiling](@ref Tiling).
        %>
        %>  \final{Tile}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:39 PM, University of Texas at Arlington<br>
        function self = Tile(template, varargin)
            [varobj, vartemp] = pm.matlab.hashmap.popKeyVal(["figure", "subplot", "template", "tiledlayout", "tileshape"], varargin);
            self = self@[pm.vis.figure.Tiling](@ref Tiling)(cell(0, 0), varobj{:});
            self.template = template;
            if ~isempty(vartemp)
                self.template.hash2comp(vartemp);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the tile to the original default settings.
        %>  Use this method when you change many attributes of the tile and
        %>  you want to clean up and go back to the default settings.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.tile.Tile.reset() % reset the tile to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:42 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)

            self.tileshape = [];
            [vartemp, varleft] = pm.matlab.hashmap.popKeyVal("template", varargin);
            if ~isempty(vartemp)
                self.hash2comp(vartemp);
                %self.template.reset();
            end
            reset@[pm.vis.figure.Tiling](@ref Tiling)(self, varleft{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the tile settings and specifications,
        %>  make the tile, and return nothing.
        %>
        %>  \details
        %>  In making the figure, this method we call the ``make()``
        %>  methods of each of the subplot objects stored in the
        %>  ``subplot`` component.
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``subplot`` component of the parent object.
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      t = pm.vis.tile.Tile.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>
        %>      t = pm.vis.tile.Tile(pm.vis.subplot.Line());
        %>      t.make()
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:44 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            self.premake(varargin{:});

            %%%% First set the tile dimensions.

            isEllipsoid = isa(self.template, "pm.vis.subplot.Ellipse") || isa(self.template, "pm.vis.subplot.Ellipse3");

            if  isEllipsoid
                nplt = size(self.template.dims, 1);
                if  size(self.template.dims, 2) ~= 2
                    help("pm.vis.tile.Ellipse3");
                    error   ( newline ...
                            + "The condition ``size(self.template.dims, 2) == 2`` must hold." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end

            else

                lencolx = 0;
                if  isprop(self.template, "colx")
                    lencolx = pm.array.len(self.template.colx);
                end

                lencoly = 0;
                if  isprop(self.template, "coly")
                    lencoly = pm.array.len(self.template.coly);
                end

                lencolz = 0;
                if  isprop(self.template, "colz")
                    lencolz = pm.array.len(self.template.colz);
                end

                lencolc = 0;
                if  isprop(self.template, "colc")
                    lencolc = pm.array.len(self.template.colc);
                end

                nplt = max([lencolx, lencoly, lencolz, lencolc]);

            end

            %%%% Define the tile shape.

            %self.tileshape = pm.matlab.hashmap.getKeyVal("tileshape", varargin);
            if ~isempty(self.tileshape)
                nrow = self.tileshape(1);
                ncol = self.tileshape(2);
                if  prod(self.tileshape) < nplt
                    help("pm.vis.tile.Tile");
                    disp("self.tileshape");
                    disp( self.tileshape );
                    disp("nplt");
                    disp( nplt );
                    error   ( newline ...
                            + "The specified tile shape must be consistent" + newline ...
                            + "with the specified data columns to visualize." + newline ...
                            + newline ...
                            );
                end
            else
                rows = nplt : -1 : 1;
                cols = ceil(nplt ./ rows);
                circ = rows + cols + (rows .* cols) - nplt;
                [~, irow] = min(circ);
                nrow = rows(irow);
                ncol = cols(irow);
            end

            %if  isempty(self.tileindex)
            %    tileindices = 1 : nplt;
            %end

            %%%% unified colorbar using tiledlayout is only possible in MARLAB R2020b and beyond.

            unicbar = ~isEllipsoid ...
                    && lencolc == 1 ...
                    && isprop(self.template, "colorbar") ...
                    && isempty(self.template.colorbar.enabled) ...
                    && "2020a" < string(version('-release'));

            %%%% Define the subplot cell matrix.
            %%%% The following code block may be improved in
            %%%% the future to avoid full data copy to subplots.

            iplt = 0;
            %tilecounter = 0;
            self.subplot = cell(nrow, ncol);
            for irow = 1 : nrow
                for icol = 1 : ncol
                    %tilecounter = tilecounter + 1;
                    %if  tilecounter == tileindices(iplt + 1)
                    iplt = iplt + 1;
                    if  nplt < iplt
                        break;
                    end
                    %self.subplot{irow, icol} = self.template;
                    copyStream = getByteStreamFromArray(self.template);
                    self.subplot{irow, icol} = getArrayFromByteStream(copyStream);
                    if  isEllipsoid
                        self.subplot{irow, icol}.dims = self.template.dims(iplt, :);
                    else
                        if  1 < lencolx
                            self.subplot{irow, icol}.colx = self.template.colx(iplt);
                        end
                        if  1 < lencoly
                            self.subplot{irow, icol}.coly = self.template.coly(iplt);
                        end
                        if  1 < lencolz
                            self.subplot{irow, icol}.colz = self.template.colz(iplt);
                        end
                        if  1 < lencolc
                            self.subplot{irow, icol}.colc = self.template.colc(iplt);
                        elseif unicbar
                            self.subplot{irow, icol}.colorbar.enabled = false;
                        elseif lencolc == 0 && isprop(self.template, "colormap") && isempty(self.template.colormap.enabled)
                            self.subplot{irow, icol}.colormap.enabled = true;
                        end
                    end
                end
                if  nplt < iplt
                    break;
                end
            end

            make@[pm.vis.figure.Tiling](@ref Tiling)(self);

            %%%% Define a single colorbar.

            if  unicbar && isprop(self.template, "colormap") && ~isempty(self.template.colormap.enabled)
                unicbar = self.template.colormap.enabled;
            end

            if  unicbar && isprop(self.template, "colorbar") && ~isempty(self.template.colorbar.enabled)
                unicbar = self.template.colorbar.enabled;
            end

            if  unicbar

                %%%% Get the start and end positions of the leftmost, lowest, rightmost, and highest axes.

                % iplt = 0;
                % positions = zeros(4, nplt);
                % for icol = 1 : ncol
                %     for irow = 1 : nrow
                %         if  pm.introspection.istype(self.subplot{irow, icol}, "pm.vis.subplot.Subplot")
                %             iplt = iplt + 1;
                %             positions(:, iplt) = self.subplot{irow, icol}.fout.axes.Position;
                %         end
                %     end
                % end
                %
                % for i = 2 : -1 : 1
                %     start(i) = min(positions(i, :));
                %     finit(i) = max(positions(i, :) + positions(i + 2, :));
                % end

                %%%% Create new invisible axes.

                kws = struct();
                for prop =  [ "colorbar" ...
                            ]
                    if  isprop(self.template, prop)
                        kws.(prop) = self.template.comp2hash(prop);
                    end
                end

                %ax = axes("position", [start, finit], "visible", "off");

                self.fout.colorbar = colorbar(kws.colorbar{:});
                self.fout.colorbar.Layout.Tile = 'east';
                dfcopy = self.template.df.copy();
                [~, colnamc] = pm.str.locname(dfcopy.Properties.VariableNames, self.template.colc);
                ylabel(self.fout.colorbar, colnamc(1));

            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the tile settings and specifications and return nothing.
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``premake()`` method.
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      f = pm.vis.tile.Tile.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{getBorder}
        %>
        %>      f = pm.vis.tile.Tile(pm.vis.Line());
        %>      f.premake("figure", {"color", "none"})
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 22 2024, 7:46 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function premake(self, varargin)

            if ~isempty(varargin)
                [varobj, vartemp] = pm.matlab.hashmap.popKeyVal(["figure", "subplot", "template", "tiledlayout", "tileshape"], varargin);
                premake@[pm.vis.figure.Tiling](@ref Tiling)(self, varobj{:});
                self.template.hash2comp(vartemp);
                %recursive = true;
                %extensible = true;
                %insensitive = true;
                %self.template = pm.matlab.hashmap.hash2comp(vartemp, self.template, insensitive, extensible, recursive);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end