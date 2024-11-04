%>  \brief
%>  This is the base class for generating objects that contain the contents of a chain
%>  file generated by a sampler of superclass [pm.sampling.Paradram](@ref Paradram).<br>
%>
%>  \details
%>  This class is meant to be primarily internally used by the ParaMonte MATLAB library samplers.<br>
%>  This class merely adds a number of visualizations to its superclass which are
%>  specific to the [pm.sampling.Paradram](@ref Paradram) sampler.<br>
%>
%>  \note
%>  See also the documentation of the class constructor [pm.sampling.FileContentsChainDRAM::FileContentsChainDRAM](@ref FileContentsChainDRAM::FileContentsChainDRAM).<br>
%>
%>  \note
%>  See also the documentation of the superclass [pm.sampling.FileContentsChainDRAM](@ref FileContentsChainDRAM).<br>
%>
%>  \note
%>  See below for information on the attributes (properties).<br>
%>
%>  \note
%>  See below for information on the class methods.<br>
%>
%>  \final{FileContentsChainDRAM}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsChainDRAM < pm.sampling.FileContentsChain

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
        %>  Return a scalar object of class [pm.sampling.FileContentsChainDRAM](@ref FileContentsChainDRAM).<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.sampling.FileContentsChainDRAM](@ref FileContentsChainDRAM).<br>
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string containing the path or web address to an external file.<br>
        %>  \param[in]  silent  :   See the corresponding argument of the superclass constructor [pm.sampling.FileContentsChain::FileContentsChain](@ref FileContentsChain::FileContentsChain).<br>
        %>                          (**optional**. The default is set by the superclass constructor [pm.sampling.FileContentsChain::FileContentsChain](@ref FileContentsChain::FileContentsChain).)
        %>  \param[in]  sep     :   See the corresponding argument of the superclass constructor [pm.sampling.FileContentsChain::FileContentsChain](@ref FileContentsChain::FileContentsChain).<br>
        %>                          (**optional**. The default is set by the superclass constructor [pm.sampling.FileContentsChain::FileContentsChain](@ref FileContentsChain::FileContentsChain).)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.FileContentsChainDRAM](@ref FileContentsChainDRAM).<br>
        %>
        %>  \interface{FileContentsChainDRAM}
        %>  \code{.m}
        %>
        %>      contents = pm.sampling.FileContentsChainDRAM(file)
        %>      contents = pm.sampling.FileContentsChainDRAM(file, [])
        %>      contents = pm.sampling.FileContentsChainDRAM(file, silent)
        %>      contents = pm.sampling.FileContentsChainDRAM(file, [], sep)
        %>      contents = pm.sampling.FileContentsChainDRAM(file, silent, sep)
        %>
        %>  \endcode
        %>
        %>  \remark
        %>  See also the documentation of the superclass the superclass constructor
        %>  [pm.sampling.FileContentsChain::FileContentsChain](@ref FileContentsChain::FileContentsChain)
        %>  for more example usage and visualizations.<br>
        %>
        %>  \example{FileContentsChainDRAM}
        %>  \include{lineno} example/sampling/FileContentsChainDRAM/main.m
        %>  \vis{FileContentsChainDRAM}
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChainDRAM/FileContentsChainDRAM.proposalAdaptation.scatter.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChainDRAM/FileContentsChainDRAM.proposalAdaptation.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChainDRAM/FileContentsChainDRAM.meanAcceptanceRate.line.png width=700
        %>  <br><br>
        %>  \image html example/sampling/FileContentsChainDRAM/FileContentsChainDRAM.burninLocation.line.png width=700
        %>
        %>  \final{FileContentsChainDRAM}
        %>
        %>  \author
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = FileContentsChainDRAM(file, silent, sep)

            if nargin < 3
                sep = [];
            end

            if nargin < 2
                silent = [];
            end

            self = self@pm.sampling.FileContentsChain(file, silent, sep);

            %%%%
            %%%% Add the ParaDRAM-specific visualizations.
            %%%%

            self.checkpoint("adding DRAM-specific visualizations to the chain object...", false);

            %%%%
            %%%% Add the visualizations for ``meanAcceptanceRate``.
            %%%%

            try

                %%%% Find the keywords in the column names.
                %%%% This avoids exact matching, which is problematic between ParaMonte 1 and 2).
                for colname = string(self.df.Properties.VariableNames)
                    if  contains(lower(colname), "mean")
                        break;
                    end
                end

                self.vis.(colname) = struct();
                self.vis.(colname).line = pm.vis.PlotLine   ( @()self.df, "coly", colname, "colc", self.slfc ...
                                                            , "ylabel", {"txt", "Mean Acceptance Rate"} ...
                                                            , "xlabel", {"txt", "Sampling Step"} ...
                                                            , "axes", {"xscale", "log"} ...
                                                            , "plot", {"linewidth", 3} ...
                                                            );

            catch me

                warning ( newline ...
                        + "Failed to create the visualizations for the " + colname + " data column of the chain file." + newline ...
                        + "Here is the error message:" + newline ...
                        + newline ...
                        + string(me.identifier) + newline + string(me.message) + newline ...
                        + newline ...
                        );

            end

            %%%%
            %%%% Add the visualizations for ``burninLocation``.
            %%%%

            try

                %%%% Find the keywords in the column names.
                %%%% This avoids exact matching, which is problematic between ParaMonte 1 and 2).
                for colname = [string(self.df.Properties.VariableNames)]
                    if  contains(lower(colname), "burnin")
                        break;
                    end
                end

                self.vis.(colname) = struct();
                self.vis.(colname).line = pm.vis.PlotLine   ( @()self.df, "coly", colname, "colc", self.slfc ...
                                                            , "ylabel", {"txt", "MCMC Burnin Location"} ...
                                                            , "xlabel", {"txt", "Sampling Step"} ...
                                                            , "axes", {"xscale", "log"} ...
                                                            , "plot", {"linewidth", 3} ...
                                                            );

            catch me

                warning ( newline ...
                        + "Failed to create the visualizations for the " + colname + " data column of the chain file." + newline ...
                        + "Here is the error message:" + newline ...
                        + newline ...
                        + string(me.identifier) + newline + string(me.message) + newline ...
                        + newline ...
                        );

            end

            %%%%
            %%%% Add the visualizations for ``proposalAdaptation``.
            %%%%

            try


                %%%% Find the keywords in the column names.
                %%%% This avoids exact matching, which is problematic between ParaMonte 1 and 2).
                for colname = [string(self.df.Properties.VariableNames)]
                    if  contains(lower(colname), "adaptation")
                        break;
                    end
                end

                self.vis.(colname) = struct();
                self.vis.(colname).line = pm.vis.PlotLine   ( @()self.df, "coly", colname, "colc", self.slfc ...
                                                            , "ylabel", {"txt", "Proposal Adaptation"} ...
                                                            , "xlabel", {"txt", "Sampling Step"} ...
                                                            , "axes", {"yscale", "log"} ...
                                                            , "plot", {"linewidth", 2} ...
                                                            );
                self.vis.(colname).scatter = pm.vis.PlotScatter ( @()self.df, "coly", colname, "colc", self.slfc ...
                                                                , "ylabel", {"txt", "Proposal Adaptation"} ...
                                                                , "xlabel", {"txt", "Sampling Step"} ...
                                                                , "axes", {"yscale", "log"} ...
                                                                );

            catch me

                warning ( newline ...
                        + "Failed to create the visualizations for the " + colname + " data column of the chain file." + newline ...
                        + "Here is the error message:" + newline ...
                        + newline ...
                        + string(me.identifier) + newline + string(me.message) + newline ...
                        + newline ...
                        );

            end

            %%%%
            %%%% Report the timing.
            %%%%

            self.checkpoint();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end