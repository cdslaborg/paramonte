%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a sample/chain file
%>  generated by a ParaMonte sampler.<br>
%>  This class is meant to be primarily internally
%>  used by the ParaMonte MATLAB library samplers.<br>
%>
%>  \note
%>  See the documentation of the class constructor
%>  [pm.sampling.FileContentsSampleDRAM::FileContentsSampleDRAM](@ref FileContentsSampleDRAM::FileContentsSampleDRAM).<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \note
%>  See also the documentation of the superclass [pm.sampling.FileContentsSample](@ref FileContentsSample).<br>
%>
%>  \final{FileContentsSampleDRAM}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 3:32 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsSampleDRAM < pm.sampling.FileContentsSample

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
        %>  Return a scalar object of class [pm.sampling.FileContentsSampleDRAM](@ref FileContentsSampleDRAM).<br>
        %>  This is the constructor of the class [pm.sampling.FileContentsSampleDRAM](@ref FileContentsSampleDRAM).<br>
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string containing the path to an external file.<br>
        %>  \param[in]  silent  :   See the corresponding argument of [pm.sampling.FileContentsRestart](@ref FileContentsRestart) class.<br>
        %>                          (**optional**. The default is set by [pm.sampling.FileContentsRestart](@ref FileContentsRestart).)
        %>  \param[in]  sep     :   The input scalar MATLAB string containing the field separator used in the file.<br>
        %>                          (**optional**, default = ``","``)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.FileContentsSampleDRAM](@ref FileContentsSampleDRAM).<br>
        %>
        %>  \interface{FileContentsSampleDRAM}
        %>  \code{.m}
        %>
        %>      contents = pm.sampling.FileContentsSampleDRAM(file)
        %>      contents = pm.sampling.FileContentsSampleDRAM(file, [])
        %>      contents = pm.sampling.FileContentsSampleDRAM(file, silent)
        %>      contents = pm.sampling.FileContentsSampleDRAM(file, [], [])
        %>      contents = pm.sampling.FileContentsSampleDRAM(file, [], sep)
        %>      contents = pm.sampling.FileContentsSampleDRAM(file, silent, [])
        %>      contents = pm.sampling.FileContentsSampleDRAM(file, silent, sep)
        %>
        %>  \endcode
        %>
        %>  \example{FileContentsSampleDRAM}
        %>  \include{lineno} example/sampling/FileContentsSampleDRAM/main.m
        %>  \vis{FileContentsSampleDRAM}
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contour.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contour.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contour.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contour3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contour3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contour3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contourf.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contourf.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.contourf.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histfit.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histfit.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histfit.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histogram.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histogram.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histogram.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histogram2.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histogram2.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.histogram2.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.line.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.line.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.line.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.line3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.line3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.line3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.lineScatter.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.lineScatter.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.lineScatter.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.lineScatter3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.lineScatter3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.lineScatter3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.scatter.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.scatter.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.scatter.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.scatter3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.scatter3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.cascade.scatter3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.contour.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.contour3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.contourf.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.histfit.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.histogram.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.histogram2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.line3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.lineScatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.lineScatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.plot.scatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.contour.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.contour3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.contourf.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.histfit.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.histogram.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.histogram2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.line3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.lineScatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.lineScatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.tile.scatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.triplex.lshc2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.triplex.lshc3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSampleDRAM/FileContentsSampleDRAM.triplex.lshcf.png width=700
        %>
        %>  \final{FileContentsSampleDRAM}
        %>
        %>  \author
        %>  \AmirShahmoradi, 7:41 PM Friday, November 8, 2024, Dallas, TX<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        function self = FileContentsSampleDRAM(file, silent, sep)

            if nargin < 3
                sep = [];
            end
            if nargin < 2
                silent = [];
            end

            self = self@pm.sampling.FileContentsSample(file, silent, sep);
            self.setstats();
            self.setvis();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end