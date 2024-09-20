%>  \brief
%>  This is the base class for generating
%>  subclass of MATLAB ``handle`` superclass whose annoying
%>  methods are forcefully hidden from the user view.<br>
%>
%>  \interface{Handle}
%>  \code{.m}
%>
%>      Handle = pm.matlab.Handle()
%>
%>  \endcode
%>
%>  \example{Handle}
%>  \include{lineno} example/matlab/Handle/main.m
%>  \output{Handle}
%>  \include{lineno} example/matlab/Handle/main.out.m
%>
%>  \final{Handle}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:31 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Handle < dynamicprops%handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>  \brief
        %>  Open the documentation of the class
        %>  of the parent object on MATLAB display.
        %>
        %>  \details
        %>  This is a dynamic method of the class [pm.matlab.Handle](@ref Handle).
        %>
        %>  \interface{doc}
        %>  \code{.m}
        %>
        %>      h = pm.matlab.Handle();
        %>      h.doc();
        %>
        %>  \endcode
        %>
        %>  \final{doc}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 11:34 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function doc(self)
            doc(class(self));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>  \brief
        %>  Print help about the class of
        %>  the parent object on MATLAB display.
        %>
        %>  \details
        %>  This is a dynamic method of the class [pm.matlab.Handle](@ref Handle).
        %>
        %>  \interface{help}
        %>  \code{.m}
        %>
        %>      self.help()
        %>
        %>  \endcode
        %>
        %>  \final{help}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 11:36 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function help(self)
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
        function setKeyVal(self, field, subfield, key, val)
            if  nargin < 4
                if ~(isfield(self, field) || isprop(self, field)) || isempty(self.(field))
                    self.(field) = subfield;
                end
            elseif  nargin < 5
                if ~(isfield(self.(field), subfield) || isprop(self.(field), subfield)) || isempty(self.(field).(subfield))
                    self.(field).(subfield) = key;
                end
            else
                if ~(isfield(self.(field).(subfield), key) || isprop(self.(field).(subfield), key)) || isempty(self.(field).(subfield).(key))
                    self.(field).(subfield).(key) = val;
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

    %>  \cond excluded

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

    %>  \endcond excluded

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end