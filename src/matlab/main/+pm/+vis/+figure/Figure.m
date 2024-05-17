classdef Figure < pm.matlab.Handle
    %
    %   This is the abstract class for generating instances of objects
    %   that contain the specifications of various types of Figures.
    %
    %   This is a generic class for generating figures containing
    %   arbitrary number of subplots (to be added by the subclasses).
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
    %       self
    %
    %           The output scalar object of class ``pm.vis.figure.Figure``.
    %
    %   Interface
    %   ---------
    %
    %       f = pm.vis.figure.Figure(varargin);
    %
    %   Attributes
    %   ----------
    %
    %       See the list of class attributes below.
    %
    %
    properties(Access = public)
        %
        %       figure
        %
        %           A MATLAB ``struct`` whose fields and their values will
        %           be passed as keyword arguments to the MATLAB intrinsic ``figure``.
        %           The following are the default components of ``figure``:
        %
        %               name
        %
        %                   Name of the figure, specified as a character vector or a string scalar.
        %
        %               color
        %
        %                   Background color, specified as an RGB triplet, a
        %                   hexadecimal color code, a color name, or a short name.
        %                   If you specify ``'none'``, the background color appears black on screen,
        %                   but if you print the figure, the background prints as though the figure
        %                   window is transparent.
        %
        %               fileName
        %
        %                   Character vector or string scalar containing the file name for
        %                   saving the figure specified as a character vector or a string scalar.
        %
        %               position
        %
        %                   Location and size of the drawable area, specified as
        %                   a vector of the form ``[left bottom width height]``.
        %                   This area excludes the figure borders, title bar, menu bar, and tool bars.
        %
        %               units
        %
        %                   Units of measurement, specified as one of the following values:
        %
        %                       'pixels' | 'normalized' | 'inches' | 'centimeters' | 'points' | 'characters'
        %
        %                   The MATLAB documentation of the intrinsic ``figure``
        %                   for the meaning of the options.
        %
        %               others
        %
        %                   See the acceptable keyword arguments of the MATLAB intrinsic ``figure()``.
        %
        %           Example usage:
        %
        %               self.figure.units = pixels;
        %               self.figure.color = "none";
        %
        %           \warning
        %
        %               Keep in mind that MATLAB keyword arguments are case-INsensitive.
        %               Hence, ensure you do not add the same keyword as multiple different fields.
        %               For example, ``figure.color`` and ``figure.Color`` are the same,
        %               and only one of the two will be processed.
        %
        figure = struct();
        %
        %   fout
        %
        %       A MATLAB ``struct`` whose fields are the outputs of
        %       various plotting tools used to make the current axis.
        %
        fout = struct();
    end

    properties(Access = protected, Hidden)
        isdryrun = true; % always false after the first make() method call or reset().
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Figure(varargin)
            self.figure = struct();
            self.reset(varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
            %
            %   Configure the figure settings and specifications,
            %   make the figure, and return nothing.
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
            %       f = pm.vis.figure.Figure.make(varargin);
            %
            %   Example
            %   -------
            %
            %       f = pm.vis.figure.Figure();
            %       f.make()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
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
            hold on;

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
            %       pm.vis.figure.Figure.reset() # reset the plot to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in premake() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.figure.alphamap = [];
            self.figure.color = [];
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

        function premake(self, varargin)
            %
            %   Configure the figure settings and specifications and return nothing.
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
            %           parent object attributes, before calling the ``premake()`` method.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       f = pm.vis.figure.Figure.premake(varargin);
            %
            %   Example
            %   -------
            %
            %       f = pm.vis.figure.Figure();
            %       f.premake("figure", {"color", "none"})
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
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

        function savefig(self, file, varargin)
            %
            %   Export the current figure to the specified external file.
            %
            %   Parameters
            %   ----------
            %
            %       file
            %
            %           The input scalar MATLAB string or char vector containing
            %           the path to the external file that will contain the visualization.
            %           The specified file extension determines its type (e.g., ``.pdf``, ``.png``).
            %           If no file extension is specified, then the default ``.png`` file extension is used.
            %
            %           (**optional**.  If ``file`` is missing or empty, first, the ``fileName`` component of
            %           the ``figure`` property of the parent object will be used as the filename if not empty.
            %           Otherwise, the figure will be exported to a file in the current working directory of
            %           MATLAB with name ``figure`` suffixed by a unique number and ``.png`` extension.)
            %
            %       varargin
            %
            %           The following optional flags are also acceptable as input string arguments.
            %           This method internally uses the venerable ``export_fig`` MATLAB package.
            %           As such, it accepts all arguments that the ``export_fig()`` function accepts.
            %           If no optional argument is present, then a default set of options determined
            %           by the ``export_fig`` library will be used.
            %
            %           The recommended set of optional flags for PNG file formats is ``"-m4 -transparent"``.
            %
            %           "-<format>"
            %
            %               The input scalar (or variadic array of) MATLAB string(s)
            %               containing the output file extension(s). Options include:
            %
            %                   '-pdf', '-eps', 'emf', '-svg', '-png', '-tif', '-jpg', '-gif', '-bmp'
            %
            %               Multiple formats can be specified, without restriction.
            %               For example: savefig('-jpg', '-pdf', '-png', ...)
            %               Note that '-tif' and '-tiff' are equivalent,
            %               and so are '-jpg' and '-jpeg'.
            %
            %           "-transparent"
            %
            %               The input option indicating that the figure background is to
            %               be made transparent (PNG, PDF, TIF, EPS, EMF formats only).
            %               Specifying this option also activates `"-noinvert"`.
            %
            %           "-nocrop"
            %
            %               The input option indicating that empty margins should not be cropped.
            %
            %           "-c[<val>,<val>,<val>,<val>]"
            %
            %               The input option indicating the crop amounts. It must be
            %               a 4-element vector of numeric values: [top, right, bottom, left]
            %               where NaN/Inf indicates auto-cropping, 0 means no cropping, any
            %               other value means cropping in pixel amounts. e.g. '-c7,15,0,NaN'
            %               Note that this option is not supported by SVG and EMF formats.
            %
            %           "-p<val>"
            %
            %               The input option to pad a border of width ``val`` to exported files,
            %               where ``val`` is either a relative size with respect to cropped image
            %               size (i.e. p=0.01 adds a 1% border). For EPS & PDF formats,
            %               ``val`` can also be integer in units of 1/72" points (abs(val)>1).
            %               ``val`` can be positive (padding) or negative (extra cropping).
            %               If used, the -nocrop flag will be ignored, i.e. the image will
            %               always be cropped and then padded. Default: 0 (i.e. no padding).
            %               Note: this option is not supported by SVG and EMF formats.
            %
            %           "-m<val>"
            %
            %               The input option ``val`` indicates the factor to magnify the figure
            %               dimensions when generating bitmap outputs (does not affect vector formats).
            %               Default: '-m1' (i.e. val=1). Note: val~=1 slows down savefig.
            %
            %           "-r<val>"
            %
            %               The input option ``val`` indicates the resolution (in pixels per inch)
            %               to export bitmap and vector outputs, without changing dimensions of
            %               the on-screen figure. Default: '-r864' (for vector output only).
            %               Note: -m option overides -r option for bitmap exports only.
            %
            %          "-native"
            %
            %               The input option indicating that the output resolution (when
            %               outputting a bitmap format) should be such that the vertical
            %               resolution of the first suitable image found in the figure is
            %               at the native resolution of that image. To specify a particular
            %               image to use, give it the tag 'export_fig_native'.
            %               Notes: This overrides any value set with the -m and -r options.
            %               It also assumes that the image is displayed front-to-parallel
            %               with the screen. The output resolution is approximate and
            %               should not be relied upon. Anti-aliasing can have adverse
            %               effects on image quality (disable with the -a1 option).
            %
            %          "-a1, -a2, -a3, -a4"
            %
            %               The input option indicating the amount of anti-aliasing (AA) to
            %               use for bitmap outputs, when GraphicsSmoothing is not available.
            %               '-a1'=no AA; '-a4'=max. Default: 3 for HG1, 1 for HG2.
            %
            %          "-<renderer>"
            %
            %               The input option to force a particular renderer (painters, opengl or
            %               [in R2014a or older] zbuffer). Default value: opengl for bitmap
            %               formats or figures with patches and/or transparent annotations;
            %               painters for vector formats without patches/transparencies.
            %
            %          "-<colorspace>"
            %
            %               The input option indicating which colorspace color figures should
            %               be saved in: RGB (default), CMYK or gray. Usage example: '-gray'.
            %               Note: CMYK is only supported in PDF, EPS and TIF formats.
            %
            %          "-q<val>"
            %
            %               The input option to vary bitmap image quality (PDF, EPS, JPG formats only).
            %               A larger val, in the range 0-100, produces higher quality and
            %               lower compression. val > 100 results in lossless compression.
            %               Default: '-q95' for JPG, ghostscript prepress default for PDF,EPS.
            %               Note: lossless compression can sometimes give a smaller file size
            %               than the default lossy compression, depending on the image type.
            %
            %          "-n<val>"
            %
            %               The input option to set minimum output image size (bitmap formats only).
            %               The output size can be specified as a single value (for both rows
            %               & cols, e.g. -n200) or comma-separated values (e.g. -n300,400).
            %               Use an Inf value to keep a dimension unchanged (e.g. -n50,inf).
            %               Use a NaN value to keep aspect ratio unchanged (e.g. -n50,nan).
            %
            %          "-x<val>"
            %
            %               The input option to set maximum output image size (bitmap formats only).
            %               The output size can be specified as a single value (for both rows
            %               & cols, e.g. -x200) or comma-separated values (e.g. -x300,400).
            %               Use an Inf value to keep a dimension unchanged (e.g. -x50,inf).
            %               Use a NaN value to keep aspect ratio unchanged (e.g. -x50,nan).
            %
            %          "-s<val>"
            %
            %               The input option to scale output image to specific size (bitmap formats only).
            %               The fixed size can be specified as a single value (for rows=cols) or
            %               comma-separated values (e.g. -s300,400). Each value can be a scalar
            %               integer (signifying pixels) or percentage (e.g. -s125%). The scaling
            %               is done last, after any other cropping/rescaling due to other params.
            %
            %          "-append"
            %
            %               The input option indicating that if the file already exists the figure is to
            %               be appended as a new page, instead of being overwritten (default).
            %               PDF, TIF & GIF output formats only (multi-image GIF = animated).
            %
            %           "-bookmark"
            %
            %               The input option to indicate that a bookmark with the name of the
            %               figure is to be created in the output file (PDF format only).
            %
            %          "-clipboard"
            %
            %               The input option to save output as an image on the system clipboard.
            %
            %          "-clipboard<:format>"
            %
            %               The input copies to clipboard in the specified format:
            %               image (default), bitmap, emf, or pdf.
            %               Notes that only -clipboard (or -clipboard:image, which is the same)
            %               applies export_fig parameters such as cropping, padding etc.
            %               -clipboard:image  create a bitmap image using export_fig processing
            %               -clipboard:bitmap create a bitmap image as-is (no auto-cropping etc.)
            %               -clipboard:emf is vector format without auto-cropping; Windows-only
            %               -clipboard:pdf is vector format without cropping; not universally supported
            %
            %          "-d<gs_option>"
            %
            %               The input option to indicate a ghostscript setting. For example,
            %               -dMaxBitmap=0 or -dNoOutputFonts (Ghostscript 9.15+).
            %
            %          "-depsc"
            %
            %               The input  option to use EPS level-3 rather than the default level-2 print
            %               device. This solves some bugs with Matlab's default -depsc2 device
            %               such as discolored subplot lines on images (vector formats only).
            %
            %          "-metadata <metaDataInfo>"
            %
            %               The input adds the specified meta-data information to the
            %               exported file (PDF format only). metaDataInfo must be either a struct
            %               or a cell array with pairs of values: {'fieldName',fieldValue, ...}.
            %               Common metadata fields: Title,Author,Creator,Producer,Subject,Keywords
            %
            %          "-nofontswap"
            %
            %               The input option to avoid font swapping. Font swapping is automatically
            %               done in vector formats (only): 11 standard Matlab fonts are
            %               replaced by the original figure fonts. This option prevents this.
            %
            %          "-font_space <char>"
            %
            %               The input option to set a spacer character for font-names that
            %               contain spaces, used by EPS/PDF. Default: ''
            %
            %          "-linecaps"
            %
            %               The input option to create rounded line-caps (vector formats only).
            %
            %          "-noinvert"
            %
            %               The input option to avoid setting figure's InvertHardcopy property to
            %               'off' during output (this solves some problems of empty outputs).
            %
            %          "-preserve_size"
            %
            %               The input option to preserve the figure's PaperSize property in output
            %               file (PDF/EPS formats only; default is to not preserve it).
            %
            %          "-options <optionsStruct>"
            %
            %               The input format-specific parameters as defined in Matlab's
            %               documentation of the imwrite function, contained in a struct under
            %               the format name. For example to specify the JPG Comment parameter,
            %               pass a struct such as this: options.JPG.Comment='abc'. Similarly,
            %               options.PNG.BitDepth=4. Only used by PNG,TIF,JPG,GIF output formats.
            %               Options can also be specified as a cell array of name-value pairs,
            %               e.g. {'BitDepth',4, 'Author','Yair'}. These options will be used
            %               by all supported output formats of the export_fig command.
            %
            %          "-silent"
            %
            %               The input option to avoid various warning and informational messages, such
            %               as version update checks, transparency or renderer issues, etc.
            %
            %          "-notify"
            %
            %               The input option to notify the user when export is done, in both a console
            %               message and a popup dialog (allow opening the exported file/folder).
            %
            %          "-regexprep <old> <new>"
            %
            %               The input replaces all occurrences of <old> (a regular expression
            %               string or array of strings; case-sensitive), with the corresponding
            %               <new> string(s), in EPS/PDF files (only). See regexp function's doc.
            %               Warning: invalid replacement can make your EPS/PDF file unreadable!
            %
            %          "-toolbar"
            %
            %               The input adds an interactive export button to the figure's toolbar
            %
            %          "-menubar"
            %
            %               The input adds an interactive export menu to the figure's menubar
            %
            %          "-contextmenu"
            %
            %               The input adds interactive export menu to figure context-menu (right-click)
            %
            %           "handle"
            %
            %               The input  handle of the figure, axes or uipanels (can be an array of handles
            %               but all the objects must be in the same figure) to be exported.
            %               Default: gcf (handle of current figure).
            %
            %           "figName"
            %
            %               The input name (title) of the figure to export (e.g. 'Figure 1' or 'My fig').
            %               Overridden by handle (if specified); Default: current figure
            %
            %   Returns
            %   -------
            %
            %       imageData
            %
            %           The output image cube of type ``uint8`` of
            %           shape ``[M, N, C]`` containing the exported figure data.
            %
            %       alpha
            %
            %           The output image matrix of shape ``[M, N]`` of alpha-matte
            %           values in the range [0,1] for the case of transparent background.
            %
            %   Interface
            %   ---------
            %
            %       f = pm.vis.figure.Figure();
            %       [imageData, alpha] = f.savefig();
            %       [imageData, alpha] = f.savefig(file);
            %       [imageData, alpha] = f.savefig(file, varargin{:});
            %
            %   Example
            %   -------
            %
            %       expoortFig(); % export the current figure with the default name.
            %       expoortFig("gridplot.pdf") % export figure to the specified PDF file.
            %       expoortFig("gridplot.png", "-m4 -transparent") % export a large png plot of magnitude 4 with transparency.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
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
                try
                    set(self.fout.figure, "color", "none");
                catch
                    failed = true;
                end
                try
                    set(gca, "color", "none");
                catch
                    failed = true;
                end
                if  failed
                    warning("Failed to set the color property of gcf to ""none"" for transparent image export.");
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
                try
                    set(self.currentFig.figure, "color", "default");
                catch
                    failed = true;
                end
                try
                    set(gca,"color","default");
                catch
                    failed = true;
                end
                if  failed
                    warning("Failed to set the color property of gca back to ""default"".");
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end