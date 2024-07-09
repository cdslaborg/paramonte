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
%>  [pm.vis.Figure](@ref Figure)<br>
%>  [pm.vis.Corner](@ref Corner)<br>
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
        %>  \image html example/vis/Cascade/Cascade.window.3.png width=700
        %>  \image html example/vis/Cascade/Cascade.window.4.png width=700
        %>
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

                nplt = max(size(self.template.subplot.rows, 1), size(self.template.subplot.colx, 1));
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

                nplt = max(1, numel(self.template.subplot.dimx), numel(self.template.subplot.dimy));

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

            self.window = cell(nplt);
            for iplt = 1 : nplt
                copyStream = getByteStreamFromArray(self.template);
                self.window{iplt} = getArrayFromByteStream(copyStream);
                if isHeatmap
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
        %>  \param[in]      files   :   The input vector of MATLAB strings or cell array of char vectors that must be
        %>                              of the same length as the length of the ``window`` component of the parent object.<br>
        %>                              containing the paths to the external files to contain the visualizations.<br>
        %>                              For more information, see the corresponding argument of the ``savefig``
        %>                              method of class [pm.vis.figure.Figure](@ref Figure).<br>
        %>                              (**optional**.  If ``files`` or any elements of it are is missing or empty,
        %>                              the default will be set by the ``savefig`` method of the corresponding
        %>                              cascade figure in the ``window`` component.)
        %>
        %>  \param[in]  varargin    :   Any optional key-val pair that is accepted by
        %>                              the ``savefig`` method of class [pm.vis.figure.Figure](@ref Figure).<br>
        %>                              For more information, see the ``savefig`` method of class [pm.vis.figure.Figure](@ref Figure).<br>
        %>
        %>  \interface{savefigs}
        %>  \code{.m}
        %>
        %>      c.savefigs();
        %>      c.savefigs(files);
        %>      c.savefigs(files, varargin);
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
        %>  \final{savefigs}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:29 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function savefigs(self, files, varargin)
            if  nargin < 2 || isempty(files)
                files = strings(length(self.window), 1);
            end
            if  length(files) ~= length(self.window)
                help("pm.vis.Cascade");
                error   ( newline ...
                        + "The condition ``length(files) == length(self.window)`` must hold." + newline ...
                        + "For more information, see the documentation displayed above." + newline ...
                        + newline ...
                        );
            end
            for iwin = 1 : length(self.window)
                self.window{iwin}.savefig(files(iwin), varargin{:});
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