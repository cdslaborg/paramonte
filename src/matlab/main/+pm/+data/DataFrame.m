%
%   This is the abstract class for generating instances of objects
%   that can contain basic attributes required for tabular read-only
%   access to a MATLAB table-compatible data stored externally.
%
%   This class is merely a convenience read-only wrapper
%   to reference external tabular data as table.
%
%   This class primarily exist to facilitate bypassing
%   the lack of references and pointers in MATLAB.
%
%       dfref
%
%           The input MATLAB 2D matrix or table containing the target dataset
%           or function handle that takes no arguments and returns the dataset.
%           Specifying a function handle is superior to specifying the dataset
%           directly, because the function handle will always allow the use of
%           the most updated version of the user table or matrix.
%           (**optional**. default = ``table(zeros(0, 0))``)
%
%>  \return
%       self
%
%           The output scalar object of class ``pm.data.DataFrame``.
%>
%>  \interface{}
%>  \code{.m}
%>  \endcode
%>
%       df = pm.data.DataFrame(dfref);
%
%   \note
%       See the class attributes descriptions below.
%
classdef DataFrame < pm.data.DataRef

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = DataFrame(dfref)
            if  nargin < 1
                dfref = [];
            end
            if  isempty(dfref)
                dfref = table(zeros(0, 0));
            end
            self = self@pm.data.DataRef(dfref);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Generate and return a table copy of the dataframe contained in the
        %   user-specified input ``dfref`` to the constructor of the parent object.
        %
        %   This class method offer the only way to access the user-specified dataframe.
        %   The underlying logic behind the use of function to access the dataframe
        %   originates from the lack of the concept of references (pointers)
        %   in the MATLAB computing language.
        %
        %       None
        %>
        %>  \return
        %>
        %       df
        %
        %           The output scalar MATLAB table a full copy of the dataframe
        %           contained in the user-specified input ``dfref`` passed
        %           to the constructor of the parent object.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       df = pm.data.DataFrame.copy()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function df = copy(self)
            df = copy@pm.data.DataRef(self);
            if ~istable(df)
                df = array2table(df);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Generate and return the number of columns in the user-specified
        %   dataframe ``dfref`` at the time of constructing the parent object.
        %
        %   This class method is a handy shorthand for ``size(self.dfref, 2)``, particularly
        %   useful for specifying a range of indices of columns in visualization tasks.
        %
        %       None
        %>
        %>  \return
        %>
        %       count
        %
        %           The output scalar MATLAB whole-number representing the number
        %           of columns in the ``dfref`` component of the parent object.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       count = pm.data.DataFrame.ncol()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function count = ncol(self)
            count = size(self.copy(), 2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Generate and return the number of rows in the user-specified
        %   dataframe ``dfref`` at the time of constructing the parent object.
        %
        %   This class method is a handy shorthand for ``size(self.dfref, 2)``,
        %   particularly useful for specifying a range of indices of rows to visualize.
        %
        %       None
        %>
        %>  \return
        %>
        %       count
        %
        %           The output scalar MATLAB whole-number representing the number
        %           of rows in the ``dfref`` component of the parent object.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       count = pm.data.DataFrame.nrow()
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function count = nrow(self)
            count = size(self.copy(), 1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %   Generate and return a natural logarithmically-spaced
        %   range of indices from the row indices of the input
        %   dataframe ``dfref`` to the parent object.
        %
        %   This method is a convenience wrapper
        %   around function ``pm.array.logrange()``.
        %
        %   \warning
        %
        %       Beware of the different order of the input arguments
        %       between this method and ``pm.array.logrange()``.
        %
        %       count
        %
        %           The input scalar MATLAB whole-number (integer)
        %           representing the maximum size of the output range.
        %           Due to rounding operation involved in creating the
        %           output range, it is impossible to prespecify the
        %           output range size, only the maximum.
        %           (**optional**, default = ``1000``)
        %
        %       start
        %
        %           The input scalar MATLAB whole-number (integer)
        %           representing the starting point of the output range.
        %           It must be a number in the range ``[1, size(self.dfref, 1)]``.
        %           Otherwise, the value ``max(1, min(start, self.nrow()))`` will be used.
        %           (**optional**, default = ``1``)
        %
        %       stop
        %
        %           The input scalar MATLAB whole-number (integer)
        %           representing the stopping point of the output range.
        %           It must be a number in the range ``[1, size(self.dfref, 1)]``.
        %           Otherwise, the value ``max(start, min(stop, self.nrow()))`` will be used.
        %           (**optional**, default = ``1``)
        %>
        %>  \return
        %>
        %       indices
        %
        %           The output vector of MATLAB real values containing
        %           the set of naturally logarithmically-spaced integer
        %           values in the specified input range.
        %
        %>  \interface{copy}
        %>  \code{.m}
        %>  \endcode
        %
        %       indices = pm.data.DataFrame()
        %       indices = pm.data.DataFrame([])
        %
        %       indices = pm.data.DataFrame([], [])
        %       indices = pm.data.DataFrame([], start)
        %       indices = pm.data.DataFrame(count, [])
        %       indices = pm.data.DataFrame(count, start)
        %
        %       indices = pm.data.DataFrame([], [], [])
        %       indices = pm.data.DataFrame(count, [], [])
        %       indices = pm.data.DataFrame([], start, [])
        %       indices = pm.data.DataFrame([], [], stop)
        %       indices = pm.data.DataFrame([], start, stop)
        %       indices = pm.data.DataFrame(count, [], stop)
        %       indices = pm.data.DataFrame(count, start, [])
        %       indices = pm.data.DataFrame(count, start, stop)
        %
        %>  \final{}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function indices = rowslog(self, count, start, stop)
            if  nargin < 4
                stop = [];
            end
            if  nargin < 3
                start = [];
            end
            if  nargin < 2
                count = [];
            end
            if  isempty(stop)
                stop = self.nrow();
            end
            if  isempty(start)
                start = 1;
            end
            if  isempty(count)
                count = 1000;
            end
            indices = pm.array.logrange(max(1, min(start, self.nrow())), max(start, min(stop, self.nrow())), count);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end