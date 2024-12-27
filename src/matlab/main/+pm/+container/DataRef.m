%>  \brief
%>  This is the abstract class for generating instances of objects
%>  that can contain basic attributes required for storing reference
%>  to a given data.<br>
%>
%>  \details
%>  This class is merely a convenience read-only
%>  wrapper to reference external data.<br>
%>
%>  \remark
%>  See the constructor documentation for example usage.<br>
%>
%>  \note
%>  This class primarily exist to facilitate bypassing
%>  the lack of references and pointers in MATLAB.<br>
%>
%>  \see
%>  [DataFrame](@ref DataFrame)<br>
%>
%>  \final
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:45 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef DataRef < pm.matlab.Handle

    properties(Access = protected, Hidden)
        %>
        %>  ``ref``
        %>
        %>  A ``protected`` and ``Hidden`` copy of the user-specified
        %>  input data (or its handle) ``ref`` to the class constructor.<br>
        %>  This is an exact copy of the user-specified function handle or data.<br>
        %>  Users are supposed to access this component only via the class
        %>  method [pm.container.DataRef.copy()](@ref DataRef::copy).<br>
        %>
        ref;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Generate an return an object of class [pm.container.DataRef](@ref DataRef)
        %>  from the input dataframe or its specified input reference.<br>
        %>
        %>  \param[in]  data    :   The input object containing the target dataset
        %>                          or function handle that takes no arguments and returns the dataset.<br>
        %>                          Specifying a function handle is superior to specifying the dataset
        %>                          directly, because the function handle will always allow the use of
        %>                          the most updated version of the user table or matrix.<br>
        %>                          (**optional**. default = ``[]``)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.container.DataRef](@ref DataRef).
        %>
        %>  \interface{DataRef}
        %>  \code{.m}
        %>
        %>      dref = pm.container.DataRef(data);
        %>
        %>  \endcode
        %>
        %>  \example{DataRef}
        %>  \include{lineno} example/container/DataRef/main.m
        %>  \output{DataRef}
        %>  \include{lineno} example/container/DataRef/main.out.m
        %>
        %>  \final{DataRef}
        %>
        %>  \author
        %>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin<br>
        function self = DataRef(data)
            if  nargin < 1
                data = [];
            end
            self.ref = data;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>
        %>  \brief
        %>  Generate and return a copy of the data (or what is returned by the function handle)
        %>  given the constructor of the parent object or the data return by the input.<br>
        %>
        %>  \details
        %>  This class method offer the only way to access the user-specified data (reference).<br>
        %>  The underlying logic behind the use of function to access the data (reference)
        %>  originates from the lack of the concept of references (pointers)
        %>  in the MATLAB computing language.<br>
        %>
        %>  \param[inout]   self    :   The **implicitly-passed** input argument representing the parent object of the method.<br>
        %>
        %>  \return
        %>  ``dfcopy``              :   The output scalar MATLAB ``table`` a full copy of the data (reference)
        %>                              contained in the user-specified input ``data`` passed
        %>                              to the constructor of the parent object.<br>
        %>
        %>  \interface{copy}
        %>  \code{.m}
        %>
        %>      dref = pm.container.DataRef(data);
        %>      dfcopy = dref.copy();
        %>
        %>  \endcode
        %>
        %>  \remark
        %>  See the constructor documentation for example usage.<br>
        %>
        %>  \final{copy}
        %>
        %>  \author
        %>  \JoshuaOsborne, Friday May 17, 2014, 6:40 PM, The University of Texas at Arlington, Dallas, TX<br>
        %>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin<br>
        function dfcopy = copy(self)
            if  isa(self.ref, 'function_handle')
                dfcopy = self.ref();
            else
                dfcopy = self.ref;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end