%>  \brief
%>  This is the abstract class for generating instances of objects
%>  that contain the specifications of various types of Figures.<br>
%>  This is a generic class for generating figures containing
%>  arbitrary number of subplots (to be added by the subclasses).
classdef Figure < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>  \brief
        %>  A MATLAB ``struct`` whose fields are the outputs of
        %>                          various plotting tools used to make the current axis.
        %>
        %>  \param[in]  figure      :   A MATLAB ``struct`` whose fields and their values will
        %>                              be passed as keyword arguments to the MATLAB intrinsic ``figure``.
        %>                              The following are the default components of ``figure``:
        %>
        %>  \param[in]  name        :   Name of the figure, specified as a character vector or a string scalar.
        %>
        %>  \param[in]  color       :   Background color, specified as an RGB triplet, a
        %>                              hexadecimal color code, a color name, or a short name.
        %>                              If you specify ``'none'``, the background color appears black on screen,
        %>                              but if you print the figure, the background prints as though the figure
        %>                              window is transparent.
        %>
        %>  \param[in]  fileName    :   Character vector or string scalar containing the file name for
        %>                              saving the figure specified as a character vector or a string scalar.
        %>
        %>  \param[in]  position    :   Location and size of the drawable area, specified as
        %>                              a vector of the form ``[left bottom width height]``.
        %>                              This area excludes the figure borders, title bar, menu bar, and tool bars.
        %>
        %>  \param[in]  units       :   Units of measurement, specified as one of the following values:<br>
        %>                                  'pixels' | 'normalized' | 'inches' | 'centimeters' | 'points' | 'characters'<br>
        %>                              The MATLAB documentation of the intrinsic ``figure``
        %>                              for the meaning of the options.
        %>
        %>  \param[in]  others      :   See the acceptable keyword arguments of the MATLAB intrinsic ``figure()``.<br>
        %>                              Example usage:<br>
        %>                                  self.figure.units = pixels;
        %>                                  self.figure.color = "none";
        %>
        %>  \warning
        %>  Keep in mind that MATLAB keyword arguments are case-INsensitive.
        %>  Hence, ensure you do not add the same keyword as multiple different fields.
        %>  For example, ``figure.color`` and ``figure.Color`` are the same,
        %>  and only one of the two will be processed.
        %>
        figure = struct();
        %>
        %>  \param[in]  fout    :   A MATLAB ``struct`` whose fields are the outputs of
        %>                          various plotting tools used to make the current axis.
        %>
        fout = struct();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = protected, Hidden)
        isdryrun = true; % always false after the first make() method call or reset().
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \return
        %>  ``self``                :   The output scalar object of class ``pm.vis.figure.Figure``.
        %>
        %>  \interface{Figure}
        %>  \code{.m}
        %>
        %>      f = pm.vis.figure.Figure(varargin);
        %>
        %>  \endcode
        %>
        %>  \note
        %>  See the list of class attributes below.
        %>
        %>  \final{Figure}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 12:06 AM, University of Texas at Arlington<br>
        function self = Figure(varargin)
            self.figure = struct();
            self.reset(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the figure settings and specifications,
        %>  make the figure, and return nothing.
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
        %>  \interface{make}
        %>  \code{.m}
        %>
        %>      f = pm.vis.figure.Figure.make(varargin);
        %>
        %>  \endcode
        %>
        %>  \example{make}
        %>
        %>      f = pm.vis.figure.Figure();
        %>      f.make()
        %>
        %>  \final{make}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 8:09 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function make(self, varargin)

            self.premake(varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: No component of ``self`` is allowed to appear to the left of assignment operator, except ``fout``.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%
            %%%% Get keyword arguments.
            %%%%

            kws = struct();
            for prop =  [ "figure" ...
                        ]
                if  isprop(self, prop)
                    kws.(prop) = self.comp2hash(prop);
                end
            end

            self.fout.figure = figure(kws.figure{:});
            %hold on; This interferes with Heatmap plots.

            %%%%
            %%%% Resize the figure to allow good default visualization.
            %%%%

            if  isempty(self.figure.outerPosition)
                maxSize = get(0, 'ScreenSize');
                figSize = self.fout.figure.OuterPosition;
                figStart = [(maxSize(3) - figSize(3)) / 2, (maxSize(4) - figSize(4)) / 2];
                set(self.fout.figure, "OuterPosition", [figStart(1:2), figSize(3:4)]);
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Reset the properties of the figure to the original default settings.
        %>  Use this method when you change many attributes of the plot and
        %>  you want to clean up and go back to the default settings.
        %>
        %>  \param[in]  varargin    :   Any ``property, value`` pair of the parent object.
        %>                              If the property is a ``struct()``, then its value must be given as a cell array,
        %>                              with consecutive elements representing the struct ``property-name, property-value`` pairs.
        %>                              Note that all of these property-value pairs can be also directly set via the
        %>                              parent object attributes, before calling the ``make()`` method.
        %>
        %>  \interface{reset}
        %>  \code{.m}
        %>
        %>      pm.vis.figure.Figure.reset() % reset the plot to the default settings.
        %>
        %>  \endcode
        %>
        %>  \final{reset}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 8:11 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function reset(self, varargin)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.figure.alphamap = [];
            self.figure.color = "white";
            self.figure.colormap = [];
            self.figure.fileName = [];
            self.figure.innerPosition = [];
            self.figure.name = [];
            self.figure.numberTitle = [];
            self.figure.outerPosition = [];
            self.figure.position = [];
            self.figure.resize = [];
            self.figure.units = [];
            self.figure.visible = [];
            self.figure.windowState = [];
            if ~isempty(varargin)
                self.hash2comp(varargin); % parse arguments
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Configure the figure settings and specifications and return nothing.
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
        %>      f = pm.vis.figure.Figure.premake(varargin);
        %>
        %>  \endcode
        %>
        %>
        %>  \example{premake}
        %>
        %>      f = pm.vis.figure.Figure();
        %>      f.premake("figure", {"color", "none"})
        %>
        %>  \final{premake}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 8:13 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function premake(self, varargin)

            if ~isempty(varargin)
                self.hash2comp(varargin); % parse arguments
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% These settings must happen here so that they can be reset every time user nullifies the values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%% figure

            %if  isfield(self.figure, "color") && pm.array.len(self.figure.color) == 0
            %    self.figure.color = "none";
            %end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function hash = comp2hash(self, comp)
            %
            %   Convert the components of the input component ``comp``
            %   of the parent object into a cell array of key-val pairs.
            %
            unique = true;
            onlyfull = true;
            excludes = {"enabled"};
            hash = pm.matlab.hashmap.struct2hash(self.(comp), excludes, unique, onlyfull);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Export the current figure to the specified external file.
        %>
        %>  \details
        %>  This method internally uses the venerable ``export_fig`` MATLAB package.<br>
        %>  As such, it accepts all arguments that the ``export_fig()`` function accepts.<br>
        %>  If no optional argument is present, then a default set of options determined
        %>  by the ``export_fig`` library will be used.<br>
        %>
        %>  \param[in]  file        :   The input scalar MATLAB string or char vector containing
        %>                              the path to the external file that will contain the visualization.
        %>                              The specified file extension determines its type (e.g., ``.pdf``, ``.png``).<br>
        %>                              If no file extension is specified, then the default ``.png`` file extension is used.<br>
        %>                              (**optional**.  If ``file`` is missing or empty, first, the ``fileName`` component of
        %>                              the ``figure`` property of the parent object will be used as the filename if not empty.<br>
        %>                              Otherwise, the figure will be exported to a file in the current working directory of
        %>                              MATLAB with name ``figure`` suffixed by a unique number and ``.png`` extension.)
        %>
        %>  \param[in]  varargin    :   The following optional flags are also acceptable as input string arguments.
        %>                              This method internally uses the venerable ``export_fig`` MATLAB package.
        %>                              As such, it accepts all arguments that the ``export_fig()`` function accepts.
        %>                              If no optional argument is present, then a default set of options determined
        %>                              by the ``export_fig`` library will be used.<br>
        %>
        %>                              The recommended set of optional flags for PNG file formats is ``"-m4 -transparent"``.<br>
        %>                              <ol>
        %>                                  <li>    ``"-<format>"``
        %>
        %>                                          The input scalar (or variadic array of) MATLAB string(s)
        %>                                          containing the output file extension(s). Options include:
        %>
        %>                                          <ol>
        %>                                              <li>    ``'-pdf'``
        %>                                              <li>    ``'-eps'``
        %>                                              <li>    ``'-emf'``
        %>                                              <li>    ``'-svg'``
        %>                                              <li>    ``'-png'``
        %>                                              <li>    ``'-tif'``
        %>                                              <li>    ``'-jpg'``
        %>                                              <li>    ``'-gif'``
        %>                                              <li>    ``'-bmp'``
        %>                                          </ol>
        %>
        %>                                          Multiple formats can be specified, without restriction.<br>
        %>                                          For example:
        %>
        %>                                          \code{.m}
        %>
        %>                                              savefig('-jpg', '-pdf', '-png')
        %>
        %>                                          \endcode
        %>
        %>                                          Note that ``'-tif'`` and ``'-tiff'`` are equivalent.<br>
        %>
        %>                                  <li>    ``"-transparent"``
        %>
        %>                                          The input option indicating that the figure background is to
        %>                                          be made transparent (PNG, PDF, TIF, EPS, EMF formats only).<br>
        %>                                          Specifying this option also activates ``"-noinvert"``.<br>
        %>
        %>                                  <li>    ``"-nocrop"``
        %>
        %>                                          The input option indicating that empty margins should not be cropped.<br>
        %>
        %>                                  <li>    ``"-c[<val>,<val>,<val>,<val>]"``
        %>
        %>                                          The input option indicating the crop amounts. It must be
        %>                                          a 4-element vector of numeric values: ``[top, right, bottom, left]``
        %>                                          where NaN/Inf indicates auto-cropping, 0 means no cropping, any
        %>                                          other value means cropping in pixel amounts. e.g. ``'-c7,15,0,NaN'``
        %>                                          Note that this option is not supported by SVG and EMF formats.<br>
        %>
        %>                                  <li>    ``"-p<val>"``
        %>
        %>                                          The input option to pad a border of width ``val`` to exported files,
        %>                                          where ``val`` is either a relative size with respect to cropped image
        %>                                          size (i.e. ``p = 0.01`` adds a ``1%`` border). For EPS & PDF formats,
        %>                                          ``val`` can also be integer in units of ``1/72`` inch points (``abs(val) > 1``).<br>
        %>                                          ``val`` can be positive (padding) or negative (extra cropping).<br>
        %>                                          If used, the ``-nocrop`` flag will be ignored, i.e. the image will
        %>                                          always be cropped and then padded. Default: ``0`` (i.e. no padding).<br>
        %>                                          Note: this option is not supported by SVG and EMF formats.<br>
        %>
        %>                                  <li>   ``"-m<val>"``
        %>
        %>                                          The input option ``val`` indicates the factor to magnify the figure
        %>                                          dimensions when generating bitmap outputs (does not affect vector formats).<br>
        %>                                          Default: '-m1' (i.e. val=1). Note: val~=1 slows down savefig.<br>
        %>
        %>                                  <li>    ``"-r<val>"``
        %>
        %>                                          The input option ``val`` indicates the resolution (in pixels per inch)
        %>                                          to export bitmap and vector outputs, without changing dimensions of
        %>                                          the on-screen figure. Default: '-r864' (for vector output only).<br>
        %>                                          Note: -m option overides -r option for bitmap exports only.<br>
        %>
        %>                                  <li>    ``"-native"``
        %>
        %>                                          The input option indicating that the output resolution (when
        %>                                          outputting a bitmap format) should be such that the vertical
        %>                                          resolution of the first suitable image found in the figure is
        %>                                          at the native resolution of that image. To specify a particular
        %>                                          image to use, give it the tag 'export_fig_native'.<br>
        %>                                          Notes: This overrides any value set with the -m and -r options.<br>
        %>                                          It also assumes that the image is displayed front-to-parallel
        %>                                          with the screen. The output resolution is approximate and
        %>                                          should not be relied upon. Anti-aliasing can have adverse
        %>                                          effects on image quality (disable with the -a1 option).<br>
        %>
        %>                                  <li>    ``"-a1, -a2, -a3, -a4"``
        %>
        %>                                          The input option indicating the amount of anti-aliasing (AA) to
        %>                                          use for bitmap outputs, when GraphicsSmoothing is not available.<br>
        %>                                          '-a1'=no AA; '-a4'=max. Default: 3 for HG1, 1 for HG2.<br>
        %>
        %>                                  <li>    ``"-<renderer>"``
        %>
        %>                                          The input option to force a particular renderer (painters, opengl or
        %>                                          [in R2014a or older] zbuffer). Default value: opengl for bitmap
        %>                                          formats or figures with patches and/or transparent annotations;
        %>                                          painters for vector formats without patches/transparencies.<br>
        %>
        %>                                  <li>    ``"-<colorspace>"``
        %>
        %>                                          The input option indicating which colorspace color figures should
        %>                                          be saved in: RGB (default), CMYK or gray. Usage example: '-gray'.<br>
        %>                                          Note: CMYK is only supported in PDF, EPS and TIF formats.<br>
        %>
        %>                                  <li>    ``"-q<val>"``
        %>
        %>                                          The input option to vary bitmap image quality (PDF, EPS, JPG formats only).<br>
        %>                                          A larger val, in the range 0-100, produces higher quality and
        %>                                          lower compression. val > 100 results in lossless compression.<br>
        %>                                          Default: '-q95' for JPG, ghostscript prepress default for PDF,EPS.<br>
        %>                                          Note: lossless compression can sometimes give a smaller file size
        %>                                          than the default lossy compression, depending on the image type.<br>
        %>
        %>                                  <li>    ``"-n<val>"``
        %>
        %>                                          The input option to set minimum output image size (bitmap formats only).<br>
        %>                                          The output size can be specified as a single value (for both rows
        %>                                          & cols, e.g. -n200) or comma-separated values (e.g. -n300,400).<br>
        %>                                          Use an Inf value to keep a dimension unchanged (e.g. -n50,inf).<br>
        %>                                          Use a NaN value to keep aspect ratio unchanged (e.g. -n50,nan).<br>
        %>
        %>                                  <li>    ``"-x<val>"``
        %>
        %>                                          The input option to set maximum output image size (bitmap formats only).<br>
        %>                                          The output size can be specified as a single value (for both rows
        %>                                          & cols, e.g. -x200) or comma-separated values (e.g. -x300,400).<br>
        %>                                          Use an Inf value to keep a dimension unchanged (e.g. -x50,inf).<br>
        %>                                          Use a NaN value to keep aspect ratio unchanged (e.g. -x50,nan).<br>
        %>
        %>                                  <li>    ``"-s<val>"``
        %>
        %>                                          The input option to scale output image to specific size (bitmap formats only).
        %>                                          The fixed size can be specified as a single value (for rows=cols) or
        %>                                          comma-separated values (e.g. -s300,400). Each value can be a scalar
        %>                                          integer (signifying pixels) or percentage (e.g. -s125%). The scaling
        %>                                          is done last, after any other cropping/rescaling due to other params.<br>
        %>
        %>                                  <li>    ``"-append"``
        %>
        %>                                          The input option indicating that if the file already exists the figure is to
        %>                                          be appended as a new page, instead of being overwritten (default).<br>
        %>                                          PDF, TIF & GIF output formats only (multi-image GIF = animated).<br>
        %>
        %>                                   <li>   ``"-bookmark"``
        %>
        %>                                          The input option to indicate that a bookmark with the name of the
        %>                                          figure is to be created in the output file (PDF format only).<br>
        %>
        %>                                  <li>    ``"-clipboard"``
        %>
        %>                                          The input option to save output as an image on the system clipboard.<br>
        %>
        %>                                  <li>    ``"-clipboard<:format>"``
        %>
        %>                                          The input copies to clipboard in the specified format:
        %>                                          image (default), bitmap, emf, or pdf.<br>
        %>                                          Notes that only ``-clipboard`` (or ``-clipboard:image``, which is the same)
        %>                                          applies export_fig parameters such as cropping, padding etc.
        %>                                          <ol>
        %>                                              <li>    ``-clipboard:image`` create a bitmap image using export_fig processing
        %>                                              <li>    ``-clipboard:bitmap`` create a bitmap image as-is (no auto-cropping etc.)
        %>                                              <li>    ``-clipboard:emf`` is vector format without auto-cropping; Windows-only.
        %>                                              <li>    ``-clipboard:pdf`` is vector format without cropping; not universally supported
        %>                                          </ol>
        %>
        %>                                  <li>    ``"-d<gs_option>"``
        %>
        %>                                          The input option to indicate a ghostscript setting. For example,
        %>                                          ``-dMaxBitmap=0`` or ``-dNoOutputFonts`` (Ghostscript 9.15+).<br>
        %>
        %>                                  <li>    ``"-depsc"``
        %>
        %>                                          The input  option to use EPS level-3 rather than the default level-2 print
        %>                                          device. This solves some bugs with MATLAB's default -depsc2 device
        %>                                          such as discolored subplot lines on images (vector formats only).<br>
        %>
        %>                                  <li>    ``"-metadata <metaDataInfo>"``
        %>
        %>                                          The input adds the specified meta-data information to the
        %>                                          exported file (PDF format only). metaDataInfo must be either a struct
        %>                                          or a cell array with pairs of values: ``{'fieldName', fieldValue, ...}``.<br>
        %>                                          Common metadata fields: Title,Author,Creator,Producer,Subject,Keywords<br>
        %>
        %>                                  <li>    ``"-nofontswap"``
        %>
        %>                                          The input option to avoid font swapping. Font swapping is automatically
        %>                                          done in vector formats (only): 11 standard MATLAB fonts are
        %>                                          replaced by the original figure fonts. This option prevents this.<br>
        %>
        %>                                  <li>    ``"-font_space <char>"``
        %>
        %>                                          The input option to set a spacer character for font-names that
        %>                                          contain spaces, used by EPS/PDF. Default: ''<br>
        %>
        %>                                  <li>    ``"-linecaps"``
        %>
        %>                                          The input option to create rounded line-caps (vector formats only).<br>
        %>
        %>                                  <li>    ``"-noinvert"``
        %>
        %>                                          The input option to avoid setting figure's InvertHardcopy property to
        %>                                          'off' during output (this solves some problems of empty outputs).<br>
        %>
        %>                                  <li>    ``"-preserve_size"``
        %>
        %>                                          The input option to preserve the figure's PaperSize property in output
        %>                                          file (PDF/EPS formats only; default is to not preserve it).<br>
        %>
        %>                                  <li>    ``"-options <optionsStruct>"``
        %>
        %>                                          The input format-specific parameters as defined in MATLAB's
        %>                                          documentation of the ``imwrite`` function, contained in a struct under
        %>                                          the format name. For example to specify the JPG Comment parameter,
        %>                                          pass a struct such as this: ``options.JPG.Comment='abc'``.<br>
        %>                                          Similarly, ``options.PNG.BitDepth = 4``.<br>
        %>                                          Only used by PNG, TIF, JPG, GIF output formats.<br>
        %>                                          Options can also be specified as a cell array of name-value pairs,
        %>                                          e.g. ``{'BitDepth', 4, 'Author', 'Yair'}``. These options will be used
        %>                                          by all supported output formats of the ``pm.vis.figure.savefig()`` function.<br>
        %>
        %>                                  <li>    ``"-silent"``
        %>
        %>                                          The input option to avoid various warning and informational messages,
        %>                                          such as version update checks, transparency or renderer issues, etc.<br>
        %>
        %>                                  <li>    ``"-notify"``
        %>
        %>                                          The input option to notify the user when export is done, in both a console
        %>                                          message and a popup dialog (allow opening the exported file/folder).<br>
        %>
        %>                                  <li>    ``"-regexprep <old> <new>"``
        %>
        %>                                          The input replaces all occurrences of <old> (a regular expression
        %>                                          string or array of strings; case-sensitive), with the corresponding
        %>                                          new string(s), in EPS/PDF files (only). See regexp function's doc.<br>
        %>                                          Warning: invalid replacement can make your EPS/PDF file unreadable!<br>
        %>
        %>                                  <li>    ``"-toolbar"``
        %>
        %>                                          The input adds an interactive export button to the figure's toolbar<br>
        %>
        %>                                  <li>    ``"-menubar"``
        %>
        %>                                          The input adds an interactive export menu to the figure's menubar<br>
        %>
        %>                                  <li>    ``"-contextmenu"``
        %>
        %>                                          The input adds interactive export menu to figure context-menu (right-click)<br>
        %>
        %>                                  <li>    ``"handle"``<br>
        %>
        %>                                          The input  handle of the figure, axes or uipanels (can be an array of handles
        %>                                          but all the objects must be in the same figure) to be exported.<br>
        %>                                          The default is ``gcf`` (handle of current figure).<br>
        %>
        %>                                  <li>    ``"figName"``
        %>
        %>                                          The input name (title) of the figure to export (e.g. ``'Figure 1'`` or ``'My fig'``).<br>
        %>                                          Overridden by handle (if specified)<br>
        %>                                          The default is the current figure.<br>
        %>                              </ol>
        %>
        %>  \return
        %>  ``imageData``       :   The output image cube of type ``uint8`` of
        %>                          shape ``[M, N, C]`` containing the exported figure data.<br>
        %>  ``alpha``           :   The output image matrix of shape ``[M, N]`` of alpha-matte
        %>                          values in the range ``[0, 1]`` for the case of transparent background.<br>
        %>
        %>  \interface{savefig}
        %>  \code{.m}
        %>
        %>      f = pm.vis.figure.Figure();
        %>      [imageData, alpha] = f.savefig();
        %>      [imageData, alpha] = f.savefig(file);
        %>      [imageData, alpha] = f.savefig(file, varargin{:});
        %>
        %>  \endcode
        %>
        %>  \example{savefig}
        %>
        %>  \code{.m}
        %>      f = pm.vis.figure.Figure();
        %>      f.savefig(); % export the current figure with the default name.
        %>      f.savefig("gridplot.pdf") % export figure to the specified PDF file.
        %>      f.savefig("gridplot.png", "-m4 -transparent") % export a large png plot of magnitude 4 with transparency.
        %>  \endcode
        %>
        %>  \final{savefig}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 8:19 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function savefig(self, file, varargin)

            if nargin < 2
                file = [];
            end
            if  pm.array.len(file) == 0
                if ~isempty(self.figure.fileName)
                    file = string(self.figure.fileName);
                else
                    if ~isempty(self.figure.name)
                        prefix = string(self.figure.name);
                    else
                        prefix = "figure";
                    end
                    fid = 0;
                    while true
                        fid = fid + 1;
                        file = prefix + "." + string(fid) + ".png";
                        if ~isfile(file)
                            break;
                        end
                    end
                end
            end

            %%%% Check the consistency of the file name.

            if ~pm.introspection.istype(file, "string", 1)
                help("pm.vis.figure.Figure");
                disp("file = ");
                disp(file);
                error   ( newline ...
                        + "The optional input argument ``file`` to method ``savefig()``" + newline ...
                        + "must be a scalar MATLAB string or char vector containing the output file name." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Focus on the target figure to export.
            %%%%

            set(0, "CurrentFigure", self.fout.figure);

            %%%%
            %%%% Ensure transparency if figure must be transparent.
            %%%%

            if  any(contains(string(varargin), "-transparent"))
                istransparent = true;
                failed = false;
                errmsg = "";
                try
                    set(self.fout.figure, "color", "none");
                catch me
                    failed = true;
                    errmsg = string(me.identifier) + " : " + string(me.message) + newline;
                end
                try
                    set(gca, "color", "none");
                catch me
                    failed = true;
                    errmsg = string(me.identifier) + " : " + string(me.message) + newline;
                end
                if  failed
                    warning ( newline + errmsg ...
                            + "Failed to set the color property of gcf to ""none"" for transparent image export." + newline ...
                            );
                end
            else
                istransparent = false;
            end

            pm.vis.figure.savefig(file, varargin{:});

            %%%%
            %%%% Revert transparency if figure ws not transparent.
            %%%%

            if  istransparent
                failed = false;
                errmsg = "";
                try
                    set(self.currentFig.figure, "color", "default");
                catch me
                    failed = true;
                    errmsg = string(me.identifier) + " : " + string(me.message) + newline;
                end
                try
                    set(gca,"color","default");
                catch me
                    failed = true;
                    errmsg = string(me.identifier) + " : " + string(me.message) + newline;
                end
                if  failed
                    warning ( newline + errmsg ...
                            + "Failed to set the color property of gca back to ""default""." + newline ...
                            );
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end