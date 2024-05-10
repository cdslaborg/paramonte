classdef (Abstract) Handle < matlab.mixin.Copyable%dynamicprops%handle
    %
    %   This is the ``Abstract`` base class for generating
    %   subclass of MATLAB ``handle`` superclass whose annoying
    %   methods are forcefully hidden from the user view.
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
    %       Handle = pm.matlab.Handle()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        function doc(self)
            %
            %   Open the documentation of the class
            %   of the parent object on MATLAB display.
            %
            %   This is a dynamic method of the class ``pm.matlab.Handle``.
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
            %       self.doc()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            doc(class(self));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        function help(self)
            %
            %   Print help about the class of
            %   the parent object on MATLAB display.
            %
            %   This is a dynamic method of the class ``pm.matlab.Handle``.
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
            %       self.help()
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            help(class(self));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)
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
                    disp("hash{i}");
                    disp( hash{i} );
                    error("The requested object property displayed above does not exist:");
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)
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
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)
        function newprop(self, prop, val)
            if ~any(strcmp(properties(self), prop))
                self.addprop(prop);
            end
            if  2 < nargin
                self.(prop) = val;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)
        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function lh = listener(varargin)
            lh = listener@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
        function delete(varargin)
            delete@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function TF = eq(varargin)
            TF = eq@handle(varargin{:});
        end
        function TF = ne(varargin)
            TF = ne@handle(varargin{:});
        end
        function TF = lt(varargin)
            TF = lt@handle(varargin{:});
        end
        function TF = le(varargin)
            TF = le@handle(varargin{:});
        end
        function TF = gt(varargin)
            TF = gt@handle(varargin{:});
        end
        function TF = ge(varargin)
            TF = ge@handle(varargin{:});
        end
        %function TF = isvalid(varargin)
        %    TF = isvalid@handle(varargin{:});
        %end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end