classdef Figure < pm.matlab.Handle
    %
    %   This is the abstract class for generating instances of objects
    %   that contain the specifications of various types of Figures.
    %
    %   Parameters
    %   ----------
    %
    %       axes
    %
    %           The input scalar MATLAB object of superclass ``pm.vis.BaseDF``.
    %
    %   Returns
    %   -------
    %
    %       self
    %
    %           The output scalar object of class ``pm.vis.figure.Base``.
    %
    %   Interface
    %   ---------
    %
    %       figure = pm.vis.figure.Base(figure);
    %
    %   Attributes
    %   ----------
    %
    %       See the list of class attributes below.
    %
    %
    properties (Access = public)
        %
        %       axes
        %
        %           A scalar object or matrix of objects of superclass ``pm.vis.axes.BaseDF``
        %           representing the (set of) axes to display in the figure.
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
        %               For example, ``figure.  color`` and ``colorbar.Color`` are the same,
        %               and only one of the two will be processed.
        %
        %   fout
        %
        %       A MATLAB ``struct`` whose fields are the outputs of
        %       various plotting tools used to make the current axis.
        %
        fout = struct();
    end

    properties (Access = protected, Hidden)
        isdryrun = true; % always false after the first make() method call or reset().
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Figure(axes)

            if  nargin < 1 || ~pm.introspection.istype(axes, "pm.vis.axes.BaseDF") || ~pm.array.len(axes)
                help("pm.vis.figure.Base");
                error   ( newline ...
                        + "The input argument ``axes`` is missing." + newline ...
                        + "For more information, see the class documentation displayed above." + newline ...
                        + newline ...
                        );
            end

            self.axes = axes;
            self.figure = struct();
            self.resetint(); % This is the subclass method!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function make(self, varargin)
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
            %       h = pm.vis.figure.Base.make(varargin);
            %
            %   Example
            %   -------
            %
            %       h = pm.vis.figure.Base(axes);
            %       h.make("xlim", [0, 1])
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

            if self.isdryrun; return; end

            %%%%%%%%%%%%%%%%
            %%%% Make plots.
            %%%%%%%%%%%%%%%%

            for i = 1 : size(self.axes, 1)
                for j = 1 : size(self.axes, 1)
                    self.axes(i, j).make()
                end
            end

        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reset(self)
            %
            %   Reset the properties of the figure to the original default settings.
            %   Use this method when you change many attributes of the plot and
            %   you want to clean up and go back to the default settings.
            %
            %   Parameters
            %   ----------
            %
            %       None
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       pm.vis.figure.Base.reset() # reset the plot to the default settings.
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self.resetint(); % call the internal reset routine.
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function resetint(self)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% RULE 0: Any non-MATLAB-default setting must be preferably set in the make() method to override user null values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for i = 1 : size(self.axes, 1)
                for j = 1 : size(self.axes, 1)
                    self.axes(i, j).reset();
                end
            end

            %%%% figure

            self.figure.color = "none";

            self.isdryrun = true;
            self.make(); % This is the subclass method!
            self.isdryrun = false;

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

        function hash2comp(self, hash)
            hashlen = length(hash);
            selfProp = string(properties(self));
            selfPropLen = length(selfProp);
            insensitive = true;
            extensible = true;
            recursive = true;
            for i = 1 : 2 : hashlen % walk through input key val.
                propertyDoesNotExist = true;
                hashItemString = string(hash{i});
                for ip = 1 : selfPropLen % walk through object prop val.
                    if  strcmpi(hashItemString, selfProp(ip))
                        propertyDoesNotExist = false;
                        if  i < hashlen % this must be here. checks for the correct pairing of key, val.
                            if  isa(hash{i + 1}, "cell") && (isstruct(self.(selfProp(ip))) || isa(self.(selfProp(ip)), "handle") || ~isempty(properties(self.(selfProp(ip)))))
                                self.(selfProp(ip)) = pm.matlab.hashmap.hash2comp(hash{i + 1}, self.(selfProp(ip)), insensitive, extensible, recursive);
                            else
                                self.(selfProp(ip)) = hash{i + 1};
                            end
                        else
                            error("The corresponding value for the property """ + string(selfProp(ip)) + """ is missing as input argument.");
                        end
                        break;
                    end
                end
                if  propertyDoesNotExist
                    error("The requested property """ + string(hash{i}) + """ does not exist.");
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setKeyVal(self, field, key, val)
            if  nargin < 4
                if  isempty(self.(field))
                    self.(field) = key;
                end
            else
                if ~isfield(self.(field), key) || isempty(self.(field).(key))
                    self.(field).(key) = val;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function newprop(self, prop, val)
            if ~any(strcmp(properties(self), prop))
                self.addprop(prop);
            end
            if  2 < nargin
                self.(prop) = val;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end