%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a progress file
%>  generated by a ParaMonte sampler.<br>
%>
%>  \details
%>  This class is meant to be primarily internally
%>  used by the ParaMonte MATLAB library samplers.<br>
%>  See the documentation of the class constructor.<br>
%>
%>  \note
%>  See below for information on the attributes (properties).<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final{FileContentsProgress}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 1:09 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsProgress < pm.io.FileContentsTabular

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
    %    %>
    %    %>  ``vis``
    %    %>  
    %    %>  The scalar MATLAB ``struct`` containing the set of
    %    %>  predefined visualizations for the output data.<br>
    %    %>
    %    vis = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Hidden)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Return a scalar object of class [pm.sampling.FileContentsProgress](@ref FileContentsProgress).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.sampling.FileContentsProgress](@ref FileContentsProgress).
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string containing the path to an external file.<br>
        %>  \param[in]  silent  :   The input scalar MATLAB logical.<br>
        %>                          if ``true``, all descriptive messages will be suppressed.<br>
        %>                          Setting this option to ``false`` is particularly useful
        %>                          in MPI-parallel simulations.<br>
        %>                          (**optional**, default = ``false``)
        %>  \param[in]  sep     :   The input scalar MATLAB string
        %>                          containing the field separator used in the file.<br>
        %>                          (**optional**, default = ``","``)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.FileContentsProgress](@ref FileContentsProgress).<br>
        %>
        %>  \interface{FileContentsProgress}
        %>  \code{.m}
        %>
        %>      contents = pm.sampling.FileContentsProgress(file)
        %>      contents = pm.sampling.FileContentsProgress(file, [])
        %>      contents = pm.sampling.FileContentsProgress(file, silent)
        %>      contents = pm.sampling.FileContentsProgress(file, [], sep)
        %>      contents = pm.sampling.FileContentsProgress(file, silent, sep)
        %>
        %>  \endcode
        %>
        %>  \final{FileContentsProgress}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 1:11 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = FileContentsProgress(file, silent, sep)
            if nargin < 3
                sep = [];
            end
            if nargin < 2
                silent = [];
            end
            self = self@pm.io.FileContentsTabular(file, silent, sep);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
