%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a chain file
%>  generated by a ParaMonte sampler.<br>
%>  This class is meant to be primarily internally
%>  used by the ParaMonte MATLAB library samplers.<br>
%>  See the documentation of the class constructor.
%>
%>  \note
%>  See below for information on the attributes (properties).
%>
%>  \note
%>  See below for information on the methods.
%>
%>  \return
%>  An object of class ``pm.sampling.FileContentsChain``.
%>
%>  \interface{FileContentsChain}
%>  \code{.m}
%>
%>      contents = pm.sampling.FileContentsChain(file)
%>      contents = pm.sampling.FileContentsChain(file, [])
%>      contents = pm.sampling.FileContentsChain(file, silent)
%>      contents = pm.sampling.FileContentsChain(file, [], sep)
%>      contents = pm.sampling.FileContentsChain(file, silent, sep)
%>
%>  \endcode
%>
%>  \final{FileContentsChain}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 1:05 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsChain < pm.sampling.FileContentsSample
    properties(Access = public)
    end

    properties(Hidden)
    end

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class ``pm.sampling.FileContentsChain``.<br>
        %>  This is the constructor of the class ``pm.sampling.FileContentsChain``.
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string
        %>                          containing the path to an external file.
        %>  
        %>  \param[in]  silent  :   See the corresponding argument of ``pm.sampling.FileContentsSample`` class.
        %>                          (**optional**. The default is set by ``pm.sampling.FileContentsSample``.)
        %>  
        %>  \param[in]  sep     :   See the corresponding argument of ``pm.sampling.FileContentsSample`` class.
        %>                          (**optional**. The default is set by ``pm.sampling.FileContentsSample``.)
        %>
        %>  \return
        %>  `self`              :   The output scalar object of class ``pm.sampling.FileContentsChain``.
        %>
        %>  \interface{FileContentsChain}
        %>  \code{.m}
        %>
        %>      contents = pm.sampling.FileContentsChain(file)
        %>      contents = pm.sampling.FileContentsChain(file, [])
        %>      contents = pm.sampling.FileContentsChain(file, silent)
        %>      contents = pm.sampling.FileContentsChain(file, [], sep)
        %>      contents = pm.sampling.FileContentsChain(file, silent, sep)
        %>
        %>  \endcode
        %>
        %>  \final{FileContentsChain}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 1:07 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = FileContentsChain(file, silent, sep)
            if nargin < 3
                sep = [];
            end
            if nargin < 2
                silent = [];
            end
            self = self@pm.sampling.FileContentsSample(file, silent, sep);
        end

    end

end