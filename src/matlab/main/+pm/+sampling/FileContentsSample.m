%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a sample/chain file
%>  generated by a ParaMonte sampler.<br>
%>  This class is meant to be primarily internally
%>  used by the ParaMonte MATLAB library samplers.<br>
%>
%>  \note
%>  See the documentation of the class constructor.<br>
%>
%>  \note
%>  See below for information on the attributes (properties).<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final{FileContentsSample}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 3:32 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsSample < pm.io.FileContentsTabular

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``ndim``    :   The scalar MATLAB integer representing the number of
        %>                  dimensions of the domain of the objective function sampled.<br>
        %>                  This integer is also the number of columns in the file that
        %>                  correspond that contain the sampled states from the domain
        %>                  of the mathematical objective function.<br>
        %>
        ndim = 0;
        %>
        %>  ``slfc``    :   The scalar MATLAB integer representing the column
        %>                  index of the dataframe component ``df`` that contains
        %>                  the natural logarithm of the objective function values
        %>                  corresponding to the sampled states next to this column,
        %>                  such that the following relationship holds.<br>
        %>                  \code{.m}
        %>                      FileContentsSample.ndim = FileContentsSample.ncol - FileContentsSample.slfc;
        %>                  \endcode
        %>                  While this column index can be readily inferred by exploring
        %>                  the contents of the dataframe component, this column index is also
        %>                  computed and explicitly offered to conveniently slice the values of
        %>                  the sampled states and their corresponding log-function values.<br>
        %>
        slfc = 0;
        %>
        %>  ``stats``   :   The scalar MATLAB object containing the set of
        %>                  computed properties of the contents of the file.<br>
        %>
        stats = [];
        %>
        %>  ``vis``     :   The scalar MATLAB ``struct`` containing the set of
        %>                  predefined visualizations for the output data.<br>
        %>
        vis = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Hidden)
        %>
        %>  ``slfcname``    :   The scalar MATLAB string representing the column
        %>                      name of the dataframe component ``df`` that contains
        %>                      the natural logarithm of the objective function values
        %>                      corresponding to the sampled states next to this column.<br>
        %>
        slfcname = "sampleLogFunc";
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %>  \brief
        %>  Return a scalar object of class [pm.sampling.FileContentsSample](@ref FileContentsSample).<br>
        %>  This is the constructor of the class [pm.sampling.FileContentsSample](@ref FileContentsSample).<br>
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string containing the path to an external file.<br>
        %>  \param[in]  silent  :   See the corresponding argument of [pm.sampling.FileContentsRestart](@ref FileContentsRestart) class.<br>
        %>                          (**optional**. The default is set by [pm.sampling.FileContentsRestart](@ref FileContentsRestart).)
        %>  \param[in]  sep     :   The input scalar MATLAB string containing the field separator used in the file.<br>
        %>                          (**optional**, default = ``","``)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.FileContentsSample](@ref FileContentsSample).<br>
        %>
        %>  \interface{FileContentsSample}
        %>  \code{.m}
        %>
        %>      contents = pm.sampling.FileContentsSample(file)
        %>      contents = pm.sampling.FileContentsSample(file, [])
        %>      contents = pm.sampling.FileContentsSample(file, silent)
        %>      contents = pm.sampling.FileContentsSample(file, [], [])
        %>      contents = pm.sampling.FileContentsSample(file, [], sep)
        %>      contents = pm.sampling.FileContentsSample(file, silent, [])
        %>      contents = pm.sampling.FileContentsSample(file, silent, sep)
        %>
        %>  \endcode
        %>
        %>  \example{FileContentsSample}
        %>  \include{lineno} example/sampling/FileContentsSample/main.m
        %>  \vis{FileContentsSample}
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contour.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contour.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contour.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contour3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contour3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contour3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contourf.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contourf.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.contourf.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histfit.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histfit.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histfit.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histogram.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histogram.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histogram.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histogram2.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histogram2.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.histogram2.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.line.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.line.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.line.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.line3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.line3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.line3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.lineScatter.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.lineScatter.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.lineScatter.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.lineScatter3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.lineScatter3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.lineScatter3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.scatter.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.scatter.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.scatter.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.scatter3.1.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.scatter3.2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.cascade.scatter3.3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.contour.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.contour3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.contourf.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.histfit.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.histogram.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.histogram2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.line3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.lineScatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.lineScatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.plot.scatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.contour.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.contour3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.contourf.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.histfit.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.histogram.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.histogram2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.line3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.lineScatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.lineScatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.tile.scatter3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.triplex.lshc2.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.triplex.lshc3.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsSample/FileContentsSample.triplex.lshcf.png width=700
        %>
        %>  \final{FileContentsSample}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 3:36 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = FileContentsSample(file, silent, sep)

            if nargin < 3
                sep = [];
            end
            if nargin < 2
                silent = [];
            end

            self = self@pm.io.FileContentsTabular(file, silent, sep);

            for icol = 1 : self.ncol
                if strcmpi(self.df.Properties.VariableNames{icol}, self.slfcname)
                    break;
                end
            end

            self.slfc = icol;
            self.ndim = self.ncol - self.slfc;
            if  self.nrow <= self.ndim
                warning ( newline ...
                        + "There are insufficient number of states in the specified file:" + newline ...
                        + newline ...
                        + pm.io.tab() + self.file + newline ...
                        + newline ...
                        + "for computing the covariance/correlation matrices. Skipping..." + newline ...
                        + newline ...
                        );
                return;
            end

            %%%%
            %%%% statistics.
            %%%%

            self.stats = struct();

            %%%% Add chain cormat.

            self.checkpoint("computing the sample correlation matrix...", false);
            self.stats.cor = pm.stats.Cor(self.df(:, self.slfc + 1 : end));
            self.checkpoint();

            %%%% Add chain covmat.

            self.checkpoint("computing the sample covariance matrix...", false);
            self.stats.cov = pm.stats.Cov(self.df(:, self.slfc + 1 : end));
            self.checkpoint();

            %%%% Add chain acf.

            self.checkpoint("computing the sample autocorrelation...", false);
            self.stats.acf = pm.stats.AutoCorr(self.df(:, self.slfc : end));
            self.checkpoint();

            self.stats.max = struct("val", [], "loc", []);
            self.stats.min = struct("val", [], "loc", []);

            %%%% The `{:,:}` slice is essential in MATLAB ~2020a.

            [self.stats.max.val, self.stats.max.loc] = max(self.df{:,:});
            [self.stats.min.val, self.stats.min.loc] = min(self.df{:,:});
            self.stats.avg = mean(self.df{:,:});
            self.stats.std = std(self.df{:,:});

            %%%%
            %%%% visualization.
            %%%%

            colf = self.slfc;
            cols = self.slfc + [0 : self.ndim];

            self.vis = struct();

            self.vis.cascade = struct();

            self.vis.cascade.line = pm.vis.CascadeLine(@()self.df, "coly", cols, "colc", colf);
            self.vis.cascade.scatter = pm.vis.CascadeScatter(@()self.df, "coly", cols, "colc", colf);
            self.vis.cascade.lineScatter = pm.vis.CascadeLineScatter(@()self.df, "coly", cols, "colc", colf);

            self.vis.cascade.line3 = pm.vis.CascadeLine3(@()self.df, "coly", cols, "colz", colf, "colc", colf);
            self.vis.cascade.scatter3 = pm.vis.CascadeScatter3(@()self.df, "coly", cols, "colz", colf, "colc", colf);
            self.vis.cascade.lineScatter3 = pm.vis.CascadeLineScatter3(@()self.df, "coly", cols, "colz", colf, "colc", colf);

            self.vis.cascade.histfit = pm.vis.CascadeHistfit(@()self.df, "colx", cols);
            self.vis.cascade.histogram = pm.vis.CascadeHistogram(@()self.df, "colx", cols);
            self.vis.cascade.histogram2 = pm.vis.CascadeHistogram2(@()self.df, "colx", cols, "coly", colf);

            self.vis.cascade.contour = pm.vis.CascadeContour(@()self.df, "colx", cols, "coly", colf);
            self.vis.cascade.contourf = pm.vis.CascadeContourf(@()self.df, "colx", cols, "coly", colf);
            self.vis.cascade.contour3 = pm.vis.CascadeContour3(@()self.df, "colx", cols, "coly", colf);

            self.vis.plot = struct();

            self.vis.plot.line = pm.vis.PlotLine(@()self.df, "coly", cols, "colc", colf);
            self.vis.plot.scatter = pm.vis.PlotScatter(@()self.df, "coly", cols, "colc", colf);
            self.vis.plot.lineScatter = pm.vis.PlotLineScatter(@()self.df, "coly", cols, "colc", colf);

            self.vis.plot.line3 = pm.vis.PlotLine3(@()self.df, "coly", cols, "colc", colf);
            self.vis.plot.scatter3 = pm.vis.PlotScatter3(@()self.df, "coly", cols, "colz", colf, "colc", colf);
            self.vis.plot.lineScatter3 = pm.vis.PlotLineScatter3(@()self.df, "coly", cols, "colz", colf, "colc", colf);

            self.vis.plot.histfit = pm.vis.PlotHistfit(@()self.df, "colx", cols);
            self.vis.plot.histogram = pm.vis.PlotHistogram(@()self.df, "colx", cols);
            self.vis.plot.histogram2 = pm.vis.PlotHistogram2(@()self.df, "colx", colf, "coly", cols);

            self.vis.plot.contour = pm.vis.PlotContour(@()self.df, "colx", cols, "coly", colf);
            self.vis.plot.contourf = pm.vis.PlotContourf(@()self.df, "colx", cols, "coly", colf);
            self.vis.plot.contour3 = pm.vis.PlotContour3(@()self.df, "colx", cols, "coly", colf);

            self.vis.tile = struct();

            self.vis.tile.line = pm.vis.TileLine(@()self.df, "coly", cols, "colc", colf);
            self.vis.tile.scatter = pm.vis.TileScatter(@()self.df, "coly", cols, "colc", colf);
            self.vis.tile.lineScatter = pm.vis.TileLineScatter(@()self.df, "coly", cols, "colc", colf);

            self.vis.tile.line3 = pm.vis.TileLine3(@()self.df, "coly", cols, "colc", colf);
            self.vis.tile.scatter3 = pm.vis.TileScatter3(@()self.df, "coly", cols, "colz", colf, "colc", colf);
            self.vis.tile.lineScatter3 = pm.vis.TileLineScatter3(@()self.df, "coly", cols, "colz", colf, "colc", colf);

            self.vis.tile.histfit = pm.vis.TileHistfit(@()self.df, "colx", cols);
            self.vis.tile.histogram = pm.vis.TileHistogram(@()self.df, "colx", cols);
            self.vis.tile.histogram2 = pm.vis.TileHistogram2(@()self.df, "colx", colf, "coly", cols);

            self.vis.tile.contour = pm.vis.TileContour(@()self.df, "colx", cols, "coly", colf);
            self.vis.tile.contourf = pm.vis.TileContourf(@()self.df, "colx", cols, "coly", colf);
            self.vis.tile.contour3 = pm.vis.TileContour3(@()self.df, "colx", cols, "coly", colf);

            self.vis.triplex.lshc2 = pm.vis.Triplex ( pm.vis.SubplotLineScatter(@()self.df, "colx", cols, "coly", cols, "colc", colf, "colorbar", {"enabled", false}) ...
                                                    , pm.vis.SubplotHistogram(@()self.df, "colx", cols) ...
                                                    , pm.vis.SubplotContour(@()self.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                                                    );

            self.vis.triplex.lshcf = pm.vis.Triplex ( pm.vis.SubplotLineScatter(@()self.df, "colx", cols, "coly", cols, "colc", colf, "colorbar", {"enabled", false}) ...
                                                    , pm.vis.SubplotHistogram(@()self.df, "colx", cols) ...
                                                    , pm.vis.SubplotContourf(@()self.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                                                    );

            self.vis.triplex.lshc3 = pm.vis.Triplex ( pm.vis.SubplotLineScatter(@()self.df, "colx", cols, "coly", cols, "colc", colf, "colorbar", {"enabled", false}) ...
                                                    , pm.vis.SubplotHistogram(@()self.df, "colx", cols) ...
                                                    , pm.vis.SubplotContour3(@()self.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                                                    );

            %if  5 < ndim
            %    self.vis.triplex.layout.tiling.position = [.1, .1, nan, nan];
            %    self.vis.triplex.layout.cbarh.position = [];
            %    self.vis.triplex.layout.cbarv.position = [];
            %    self.vis.triplex.layout.tiling.tile.width = [];
            %end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % %>  \brief
        % %>  Visualize, quickly and minimally, a selected subset of the contents of the parent object and return nothing.<br>
        % %>
        % %>  \details
        % %>  All visualizations created belong to the ``vis`` component of the parent object.<br>
        % %>  This method is merely convenience to automate some of the most desired visualizations in the ``vis`` component.<br>
        % %>
        % %>  \param[in]  self    :   The input/output parent object of class [pm.sampling.FileContentsSample](@ref FileContentsSample)
        % %>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
        % %>
        % %>  \interface{visit}
        % %>  \code{.m}
        % %>
        % %>      sample = pm.sampling.FileContentsSample(file);
        % %>      sample.visit();
        % %>
        % %>  \endcode
        % %>
        % %>  \warning
        % %>  This method is to be only used for post-processing of the output
        % %>  chain file(s) of an already finished simulation. It is NOT meant to
        % %>  be called by all processes in parallel mode, although it is possible.<br>
        % %>
        % %>  \example{visit}
        % %>  \include{lineno} example/sampling/Paradram/visit/main.m
        % %>  \include{lineno} example/sampling/Paradram/visit/main.out.m
        % %>  \vis{visit}
        % %>  \image html example/sampling/FileContentsSample/visit/FileContentsSample.visit.domain.png width=700
        % %>  <br><br>
        % %>  \image html example/sampling/FileContentsSample/visit/FileContentsSample.visit.traceplot.png width=700
        % %>  <br><br>
        % %>  \image html example/sampling/FileContentsSample/visit/FileContentsSample.visit.proposalAdaptation.png width=700
        % %>
        % %>  \final{visit}
        % %>
        % %>  \author
        % %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        % function visit(self, pattern, sep)
        % 
        % end
        % 
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
