%>  \brief
%>  This is the abstract class for generating instances of objects
%>  that contain the specifications of a cascade of plots.<br>
%>  This is a generic class for generating multiple separate
%>  cascading figures each containing a single plot (axes).
classdef Cascade < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  \param  window      :   The scalar object of superclass ``pm.vis.cascade.Cascade``
        %>                          representing the cascade of plots to display.
        %>
        window = [];
        %>
        %>  \param  template    :   The scalar object of superclass ``pm.vis.cascade.Cascade``
        %>                          representing the template of the cascade of plots to display.
        %>
        template = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %>
        %>  \param[in]  template    :   The input scalar object of superclass ``pm.vis.cascade.Cascade``.
        %>                              It serves as the template based upon which
        %>                              the cascade of plots are constructed.
        %>  
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>  
        %>  \note
        %>  The input ``varargin`` can also contain the components
        %>  of the ``plot`` component of the parent object.
        %>
        %>  \return
        %>  `self`                  :   The output scalar object of class ``pm.vis.cascade.Cascade``.
        %>
        %>  \interface{Cascade}
        %>  \code{.m}
        %>
        %>      plot = pm.vis.cascade.Cascade(template);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the list of class attributes below,
        %>  also those of the superclass ``pm.vis.figure.Handle``.
        %>  
        %>  \final{Cascade}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:20 AM, University of Texas at Arlington<br>

        function self = Cascade(template, varargin)
            self.template = template;
            self.reset(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the cascade of plots to the original default settings.
        %>  Use this method when you change many attributes of the plot and
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
        %>  of the ``template`` component of the parent object.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.cascade.Cascade.reset() # reset the plot to the default settings.
        %>
        %>  \endcode
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:25 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
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
        %>  make the cascade of plots, and return nothing.
        %>
        %>  \details
        %>  In making the figure, this method we call the ``make()``
        %>  methods of each of the plot objects stored in the
        %>  ``window`` component.
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
        %>  of the ``template`` component of the parent object.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      p = pm.vis.cascade.Cascade.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>
        %>      p = pm.vis.cascade.Cascade("coly", 1:3);
        %>      p.make()
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:27 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            self.premake(varargin{:});

            %%%% First set the cascade dimensions.

            isEllipsoid = isa(self.template, "pm.vis.plot.Ellipse") || isa(self.template, "pm.vis.plot.Ellipse3");

            if  isEllipsoid
                nplt = size(self.template.subplot.dims, 1);
                if  size(self.template.subplot.dims, 2) ~= 2
                    help("pm.vis.plot.Ellipse3");
                    error   ( newline ...
                            + "The condition ``size(self.template.dims, 2) == 2`` must hold." + newline ...
                            + "For more information, see the documentation displayed above." + newline ...
                            + newline ...
                            );
                end

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
                if  isEllipsoid
                    self.window{iplt}.subplot.dims = self.template.subplot.dims(iplt, :);
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
        %>  Export the current cascade of figures to the specified external files.
        %>
        %>  \note
        %>  This method is merely a wrapper around
        %>  all ``savefig`` methods of the figure objects in
        %>  the ``window`` component of the parent cascade object.
        %>
        %>  \param[in]  files   :   The input vector of MATLAB strings or cell array of char vectors that must be
        %>                          of the same length as the length of the ``window`` component of the parent object.
        %>                          containing the paths to the external files to contain the visualizations.
        %>                          For more information, see the corresponding argument of the ``savefig``
        %>                          method of class ``pm.vis.figure.Figure``.
        %>                          (**optional**.  If ``files`` or any elements of it are is missing or empty,
        %>                          the default will be set by the ``savefig`` method of the corresponding
        %>                          cascade figure in the ``window`` component.)
        %>  
        %>  \param[in]  varargin    :   Any optional key-val pair that is accepted by
        %>                              the ``savefig`` method of class ``pm.vis.figure.Figure``.
        %>                              For more information, see the ``savefig`` method of class ``pm.vis.figure.Figure``.
        %>
        %>  \return
        %>  `None`
        %>
        %>  \interface{savefigs}
        %>  \code{.m}
        %>
        %>      c.savefigs();
        %>      c.savefigs(files);
        %>      c.savefigs(files, varargin{:});
        %>  \endcode
        %>
        %>
        %>  \example{getBorder}
        %>
        %>      c = pm.vis.cascade.Cascade(pm.vis.plot.Line());
        %>      c.savefigs(); % export the current figure with the default name.
        %>      c.savefigs("gridplot.pdf") % export figure to the specified PDF file.
        %>      c.savefigs("gridplot.png", "-m4 -transparent") % export a large png plot of magnitude 4 with transparency.
        %>
        %>  \final{savefigs}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:29 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function savefigs(self, files, varargin)
            if  nargin < 2 || isempty(files)
                files = strings(length(self.window), 1);
            end
            if  length(files) ~= length(self.window)
                help("pm.vis.cascade.Cascade");
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
        %>  Configure the cascade of plots settings and specifications and return nothing.
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
        %>  \return
        %>  `None`
        %>
        %>  \interface{premake}
        %>  \code{.m}
        %>
        %>      c = pm.vis.cascade.Cascade.premake(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{premake}
        %>
        %>      c = pm.vis.cascade.Cascade();
        %>      c.premake("figure", {"color", "none"})
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:31 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
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