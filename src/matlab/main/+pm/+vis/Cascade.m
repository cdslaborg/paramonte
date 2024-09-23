%>  \brief
%>  This is the abstract class for generating instances of objects
%>  that contain the specifications of a cascade of figures (plots).<br>
%>
%>  \details
%>  This is a generic class for generating multiple separate
%>  cascading figures each containing a single plot (axes).<br>
%>
%>  \devnote
%>  While dynamic addition of class attributes is not ideal, the current
%>  design was deemed unavoidable and best, given the constraints of the
%>  MATLAB language and visualization tools.<br>
%>
%>  \note
%>  See the list of class attributes below,
%>  also those of the superclass [pm.matlab.Handle](@ref Handle).<br>
%>
%>  \see
%>  [pm.vis.Cascade](@ref Cascade)<br>
%>  [pm.vis.Subplot](@ref Subplot)<br>
%>  [pm.vis.Triplex](@ref Triplex)<br>
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Plot](@ref Plot)<br>
%>  [pm.vis.Tile](@ref Tile)<br>
%>
%>  \final
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Cascade < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``window``
        %>
        %>  The cell array of scalar objects of superclass [pm.vis.Plot](@ref Plot)
        %>  representing the cascade of plots to display.<br>
        %>
        window = [];
        %>
        %>  ``template``
        %>
        %>  The scalar object of superclass [pm.vis.Plot](@ref Plot)
        %>  representing the template of the cascade of plots to display.<br>
        %>
        template = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Construct and return an object of class [pm.vis.Cascade](@ref Cascade).<br>
        %>
        %>  \details
        %>  This is the custom constructor of the class [pm.vis.Cascade](@ref Cascade).<br>
        %>
        %>  \param[in]  template    :   The input scalar object of superclass [pm.vis.Plot](@ref Plot).<br>
        %>                              It serves as the template based upon which the cascade of plots are constructed.<br>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the [pm.vis.Cascade.make()](@ref Cascade::make) method.<br>
        %>
        %>  \return
        %>  ``self``                :   The scalar output object of class [pm.vis.Cascade](@ref Cascade).<br>
        %>
        %>  \interface{Cascade}
        %>  \code{.m}
        %>
        %>      plot = pm.vis.Cascade(template);
        %>      plot = pm.vis.Cascade(template, varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``plot`` component of the parent object.<br>
        %>
        %>  \example{Cascade}
        %>  \include{lineno} example/vis/Cascade/main.m
        %>  \vis{Cascade}
        %>  \image html example/vis/Cascade/Cascade.window.1.png width=700
        %>  \image html example/vis/Cascade/Cascade.window.2.png width=700
        %>  \image html example/vis/Cascade/Cascade.window.3.png width=700
        %>
        %>  \final{Cascade}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:20 AM, University of Texas at Arlington<br>
        %>  \AmirShahmoradi, July 5 2024, 1:07 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function self = Cascade(template, varargin)
            self.template = template;
            self.reset(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the cascade of plots to the original default settings.<br>
        %>
        %>  \details
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.Cascade](@ref Cascade)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the [pm.vis.Cascade.make()](@ref Cascade::make) method.<br>
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``template`` component of the parent object.
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.Cascade.reset() % reset the plot to the default settings.
        %>      pm.vis.Cascade.reset(varargin) % reset the plot to the default settings and to those specified via ``varargin``.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:25 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)

            if ~isempty(varargin)
                [varpar, varsub] = pm.matlab.hashmap.popKeyVal(["template", "window"], varargin);
                if ~isempty(varpar)
                    self.hash2comp(varpar);
                end
                if ~isempty(varsub)
                    self.template.reset(varsub{:});
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the cascade settings and specifications,
        %>  make the cascade of plots, and return nothing.<br>
        %>
        %>  \details
        %>  In making the figure, this method calls the ``make()``
        %>  methods of each of the plot objects stored in the
        %>  ``window`` component.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.Cascade](@ref Cascade)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the parent object attributes,
        %>                                  before calling the [pm.vis.Cascade.make()](@ref Cascade::make) method.<br>
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      p = pm.vis.Cascade.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``template`` component of the parent object.<br>
        %>
        %>  \example{Cascade}
        %>  \include{lineno} example/vis/Cascade/main.m
        %>  \vis{Cascade}
        %>  \image html example/vis/Cascade/Cascade.window.1.png width=700
        %>  \image html example/vis/Cascade/Cascade.window.2.png width=700
        %>
        %>  \image html example/vis/Cascade/Cascade.window.3.png width=700
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:27 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            self.premake(varargin{:});

            %%%% First set the cascade dimensions.

            isHeatmap = isa(self.template, "pm.vis.PlotHeatmap");
            isEllipsoid = isa(self.template, "pm.vis.PlotEllipse") || isa(self.template, "pm.vis.PlotEllipse3");

            if  isHeatmap

                %%%% Set the data for the Heatmap.

                hdf = struct();
                hdf.colx = self.template.subplot.colx;
                hdf.rows = self.template.subplot.rows;
                if ~iscell(self.template.subplot.colx) && ~isempty(self.template.subplot.colx)
                    hdf.colx = {hdf.colx};
                end
                if ~iscell(self.template.subplot.rows) && ~isempty(self.template.subplot.rows)
                    hdf.rows = {hdf.rows};
                end

                %%%% Test the data size consistency.

                nplt = max(numel(hdf.colx), numel(hdf.rows));
                if (nplt ~= size(self.template.subplot.rows, 1) && size(self.template.subplot.rows, 1) > 1) ...
                || (nplt ~= size(self.template.subplot.colx, 1) && size(self.template.subplot.colx, 1) > 1)
                    help("pm.vis.CascadeHeatmap");
                    disp("size(self.template.subplot.rows)");
                    disp( size(self.template.subplot.rows) );
                    disp("size(self.template.subplot.colx)");
                    disp( size(self.template.subplot.colx) );
                    error   ( newline ...
                            + "The size of the ``rows`` and ``colx`` attributes of the ``subplot`` component" + newline ...
                            + "of the input template ``Plot`` object must equal along the first dimension." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end

            elseif isEllipsoid

                nplt = max([1, numel(self.template.subplot.dimx), numel(self.template.subplot.dimy)]);

            else

                lencolx = 0;
                if  isprop(self.template.subplot, "colx")
                    lencolx = pm.array.len(self.template.subplot.colx);
                end

                lencoly = 0;
                if  isprop(self.template.subplot, "coly")
                    lencoly = pm.array.len(self.template.subplot.coly);
                end

                lencolz = 0;
                if  isprop(self.template.subplot, "colz")
                    lencolz = pm.array.len(self.template.subplot.colz);
                end

                lencolc = 0;
                if  isprop(self.template.subplot, "colc")
                    lencolc = pm.array.len(self.template.subplot.colc);
                end

                nplt = max([lencolx, lencoly, lencolz, lencolc]);

            end

            %%%% Define the cascade.
            %%%% The following code block may be improved in
            %%%% the future to avoid full data copy to subplots.

            self.window = cell(nplt, 1);
            for iplt = 1 : nplt
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% The byte stream approach leads to serious problems with
                %%%% figures when generated from within sampler components.
                %byteStream = getByteStreamFromArray(self.template);
                %self.window{iplt} = getArrayFromByteStream(byteStream);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                self.window{iplt} = pm.matlab.copy(self.template, eval(string(class(self.template)) + "()"));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if  isHeatmap
                    if ~isempty(hdf.colx)
                        self.window{iplt}.subplot.colx = hdf.colx{min(iplt, numel(hdf.colx))};
                    end
                    if ~isempty(hdf.rows)
                        self.window{iplt}.subplot.rows = hdf.rows{min(iplt, numel(hdf.rows))};
                    end
                elseif isEllipsoid
                    if ~isempty(self.template.subplot.dimx)
                        self.window{iplt}.subplot.dimx = self.template.subplot.dimx(min(iplt, numel(self.template.subplot.dimx)));
                    end
                    if ~isempty(self.template.subplot.dimy)
                        self.window{iplt}.subplot.dimy = self.template.subplot.dimy(min(iplt, numel(self.template.subplot.dimy)));
                    end
                else
                    if  1 < lencolx
                        self.window{iplt}.subplot.colx = self.template.subplot.colx(iplt);
                    end
                    if  1 < lencoly
                        self.window{iplt}.subplot.coly = self.template.subplot.coly(iplt);
                    end
                    if  1 < lencolz
                        self.window{iplt}.subplot.colz = self.template.subplot.colz(iplt);
                    end
                    if  1 < lencolc
                        self.window{iplt}.subplot.colc = self.template.subplot.colc(iplt);
                    end
                end
                self.window{iplt}.make();
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Export the current cascade of figures to the specified external files.<br>
        %>
        %>  \details
        %>  This method is merely a wrapper around
        %>  all ``savefig`` methods of the figure objects in
        %>  the ``window`` component of the parent cascade object.<br>
        %>
        %>  \param[inout]   self    :   The input/output parent object of class [pm.vis.Cascade](@ref Cascade)
        %>                              which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      file    :   The input **vector** of MATLAB strings or cell array of char vectors
        %>                              containing the paths to the external output figure files.<br>
        %>                              It must be of the same length as the ``window`` component of the parent object.<br>
        %>                              For more information, see the corresponding argument ``file`` of the ``savefig``
        %>                              method of class [pm.vis.figure.Figure](@ref Figure).<br>
        %>                              (**optional**.  If ``file`` or any elements of it is missing or empty,
        %>                              the default will be set by the ``savefig`` method of the corresponding
        %>                              cascade figure in the ``window`` component.)
        %>
        %>  \param[in]  varargin    :   Any optional key-val pair that is accepted by
        %>                              the ``savefig`` method of class [pm.vis.figure.Figure](@ref Figure).<br>
        %>                              For more information, see the ``savefig`` method of class [pm.vis.figure.Figure](@ref Figure).<br>
        %>
        %>  \interface{savefig}
        %>  \code{.m}
        %>
        %>      c.savefig();
        %>      c.savefig(file);
        %>      c.savefig(file, varargin);
        %>
        %>  \endcode
        %>
        %>  \example{Cascade}
        %>  \include{lineno} example/vis/Cascade/main.m
        %>  \vis{Cascade}
        %>  \image html example/vis/Cascade/Cascade.window.1.png width=700
        %>  \image html example/vis/Cascade/Cascade.window.2.png width=700
        %>  \image html example/vis/Cascade/Cascade.window.3.png width=700
        %>
        %>  \final{savefig}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:29 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function savefig(self, file, varargin)
            % deleted doc for ``file`` behavior:
            %                              It must be either,
            %                              <ol>
            %                                  <li>    of the same length as the ``window`` component of the parent object, or,<br>
            %                                  <li>    scalar string or vector of length one.<br>
            %                                          If so, the specified single value for ``file`` will be used a prefix
            %                                          for all output figures and each figure will be suffixed with ``.i.png``
            %                                          where ``i`` is replaced with the figure number starting from ``1``.<br>
            %                                          **Beware that any existing figure file with the same name will be overwritten.**<br>
            %                              </ol>
            if  nargin < 2 || isempty(file)
                file = strings(numel(self.window), 1);
            end
            if  numel(file) ~= numel(self.window)
                help("pm.vis.Cascade");
                disp("numel(self.window)");
                disp( numel(self.window) );
                disp("numel(file)");
                disp( numel(file) );
                error   ( newline ...
                        + "The condition ``numel(file) == numel(self.window)`` must hold." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            %%%%
            %%%%    Starting from the last is important
            %%%%    to ensure all figures have been already created,
            %%%%    otherwise, the next figure to be generate will
            %%%%    become active over the previous images that
            %%%%    are being saved, corrupting the export.
            %%%%
            %%%%    \todo
            %%%%    Assuming all figures are generated sequentially from the first to the last,
            %%%%    the following fix to export figures from the last to the first should work.
            %%%%    However, this is only a temporary fix. An ideal solution must ensure all
            %%%%    figures have been generated before starting to export any figure.
            %%%%
            for iwin = numel(self.window) : -1 : 1
                self.window{iwin}.savefig(file(iwin), varargin{:});
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the cascade of plots settings and specifications and return nothing.<br>
        %>
        %>  \warning
        %>  This method has side-effects by manipulating
        %>  the existing attributes of the parent object.<br>
        %>
        %>  \param[inout]   self        :   The input/output parent object of class [pm.vis.Cascade](@ref Cascade)
        %>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]      varargin    :   Any ``property, value`` pair of the parent object.<br>
        %>                                  If the property is a ``struct()``, then its value must be given as a cell array,
        %>                                  with consecutive elements representing the struct ``property-name, property-value`` pairs.<br>
        %>                                  Note that all of these property-value pairs can be also directly set via the
        %>                                  parent object attributes, before calling the ``premake()`` method.<br>
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      c = pm.vis.Cascade.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{premake}
        %>  \code{.m}
        %>
        %>      cv = pm.vis.Cascade(pm.vis.PlotLine(rand(100, 10)));
        %>      c.premake("figure", {"color", "white"})
        %>      cv.make("coly", 1:3)
        %>
        %>  \endcode
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:31 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function premake(self, varargin)

            if ~isempty(varargin)
                [varpar, varsub] = pm.matlab.hashmap.popKeyVal(["template", "window"], varargin);
                if ~isempty(varpar)
                    self.hash2comp(varpar);
                end
                if ~isempty(varsub)
                    self.template.premake(varsub{:});
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end