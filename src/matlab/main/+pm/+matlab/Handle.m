classdef (Abstract) Handle < dynamicprops%handle
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

    methods(Hidden)
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

end