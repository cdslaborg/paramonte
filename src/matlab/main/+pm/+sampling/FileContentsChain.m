%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a chain file
%>  generated by a ParaMonte sampler.<br>
%>
%>  \details
%>  This class is meant to be primarily internally
%>  used by the ParaMonte MATLAB library samplers.<br>
%>  See the documentation of the class constructor.<br>
%>
%>  \note
%>  See below for information on the attributes (properties).
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final{FileContentsChain}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 1:05 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsChain < pm.sampling.FileContentsSample

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Hidden)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Return a scalar object of class [pm.sampling.FileContentsChain](@ref FileContentsChain).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.sampling.FileContentsChain](@ref FileContentsChain).<br>
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string containing the path or web address to an external file.<br>
        %>  \param[in]  silent  :   See the corresponding argument of [pm.sampling.FileContentsSample](@ref FileContentsSample) class.<br>
        %>                          (**optional**. The default is set by [pm.sampling.FileContentsSample](@ref FileContentsSample).)
        %>  \param[in]  sep     :   See the corresponding argument of [pm.sampling.FileContentsSample](@ref FileContentsSample) class.<br>
        %>                          (**optional**. The default is set by [pm.sampling.FileContentsSample](@ref FileContentsSample).)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.FileContentsChain](@ref FileContentsChain).<br>
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
        %>  \example{FileContentsChain}
        %>  \include{lineno} example/sampling/FileContentsChain/main.m
        %>  \vis{FileContentsChain}
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contour.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contour.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contour.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contour3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contour3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contour3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contourf.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contourf.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.contourf.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histfit.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histfit.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histfit.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histogram.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histogram.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histogram.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histogram2.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histogram2.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.histogram2.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.line.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.line.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.line.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.line3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.line3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.line3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.lineScatter.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.lineScatter.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.lineScatter.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.lineScatter3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.lineScatter3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.lineScatter3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.scatter.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.scatter.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.scatter.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.scatter3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.scatter3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.cascade.scatter3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.contour.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.contour3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.contourf.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.histfit.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.histogram.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.histogram2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.line3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.lineScatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.lineScatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.plot.scatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.contour.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.contour3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.contourf.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.histfit.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.histogram.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.histogram2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.line3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.lineScatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.lineScatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.tile.scatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.triplex.lshc2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.triplex.lshc3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChain/FileContentsChain.triplex.lshcf.png width=700
        %>
        %>  \final{FileContentsChain}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 1:07 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end