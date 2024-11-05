%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a report file
%>  generated by a ParaMonte sampler.<br>
%>
%>  \details
%>  This class is meant to be primarily internally
%>  used by the ParaMonte MATLAB library samplers.<br>
%>  See the documentation of the class constructor.
%>
%>  \note
%>  See below for information on the attributes (properties).<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final{FileContentsReport}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 1:13 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContentsReport < pm.io.FileContents

    properties(Access = public)
        %>
        %>  ``stats``       :   The scalar MATLAB ``struct`` containing the set of
        %>                      computed properties extracted from the report file.<br>
        %>
        stats = struct();
        %>
        %>  ``contents``    :   The scalar MATLAB string containing the entire
        %>                      contents of the report file with all Carriage Return
        %>                      characters removed (relevant only to Windows OS).<br>
        %>
        contents = [];
        %>
        %>  ``lineList``    :   The vector of MATLAB strings containing the set of
        %>                      all lines in the report file with all Carriage Return
        %>                      and New Line characters removed.<br>
        %>
        lineList = [];
        %>
        %>  ``vis``         :   The scalar MATLAB ``struct`` containing the set of
        %>                      predefined visualizations for the output data.<br>
        %>
        vis = [];
        %>
        %>  ``banner``      :   The scalar MATLAB ``string`` containing the ParaMonte
        %>                      library banner as appearing in the report file.<br>
        %>
        banner = "";
        % %
        % %   setup
        % %
        % %       The scalar MATLAB ``struct`` containing the sampler
        % %       setup information extracted from the report file.
        % %
        % setup = struct();
        % %
        % %   spec
        % %
        % %       The scalar MATLAB ``struct`` containing the set of
        % %       simulation specifications extracted from the report file.
        % %
        % spec = struct();
    end

    properties(Hidden)
        %>
        %>  ``lineListLen``
        %>
        %>  This is an internal class variable inaccessible to the end users.<br>
        %>
        lineListLen = [];
        %>
        %>  ``indentLen``
        %>
        %>  The scalar MATLAB integer representing the number of indentation
        %>  characters at the beginning of each description line in the report file.<br>
        %>  This is an internal class variable inaccessible to the end users.<br>
        %>
        indentLen = 4; % indent length of the records
        %>
        %>  ``dsymLen``
        %>
        %>  The scalar MATLAB integer representing the minimum length of two of decoration symbols.<br>
        %>  This is an internal class variable inaccessible to the end users.<br>
        %>
        dsymLen = 2;
        %>
        %>  ``dsym``
        %>
        %>  The scalar MATLAB string representing the decoration symbol used in
        %>  the report file, to be determined at runtime (currently ``%``).<br>
        %>  This is an internal class variable inaccessible to the end users.<br>
        %>
        dsym = '';
        %>
        %>  ``prefix``
        %>
        %>  The scalar MATLAB string representing the prefix used in the description lines of the report file.<br>
        %>  This is an internal class variable inaccessible to the end users.<br>
        %>
        prefix = ' - NOTE: ';
        %>
        %>  ``method``
        %>
        %>  The scalar MATLAB string representing the sample name.<br>
        %>  This is an internal class variable inaccessible to the end users.<br>
        %>
        method = '';
        %iend = 0;
    end

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class [pm.sampling.FileContentsReport](@ref FileContentsReport).<br>
        %>  \details
        %>  This is the constructor of the class [pm.sampling.FileContentsReport](@ref FileContentsReport).<br>
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string containing the path to an external report file.<br>
        %>  \param[in]  silent  :   See the corresponding argument of [pm.io.FileContents](@ref FileContents) class.<br>
        %>                          (**optional**. The default is set by [pm.io.FileContents](@ref FileContents).)
        %>  \param[in]  method  :   The input scalar MATLAB string
        %>                          containing the sampling method name.<br>
        %>                          The input value must be any of the following:<br>
        %>                          <ol>
        %>                              <li>    ``"ParaDRAM"``
        %>                              <li>    ``"ParaDISE"``
        %>                              <li>    ``"ParaNest"``
        %>                          </ol>
        %>                          (**optional**. If missing, some of the report file contents may not be properly parsed.)
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.FileContentsReport](@ref FileContentsReport).<br>
        %>
        %>  \interface{FileContentsReport}
        %>  \code{.m}
        %>
        %>      contents = pm.sampling.FileContentsReport(file)
        %>      contents = pm.sampling.FileContentsReport(file, [])
        %>      contents = pm.sampling.FileContentsReport(file, silent)
        %>      contents = pm.sampling.FileContentsReport(file, [], [])
        %>      contents = pm.sampling.FileContentsReport(file, silent, [])
        %>      contents = pm.sampling.FileContentsReport(file, silent, method)
        %>
        %>  \endcode
        %>
        %>  \final{FileContentsReport}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 1:30 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = FileContentsReport(file, silent, method)

            if  nargin < 2
                silent = [];
            end
            self = self@pm.io.FileContents(file, silent);
            if  2 < nargin
                self.method = convertStringsToChars(method);
                self.prefix = [self.method, self.prefix];
            end

            %%%%
            %%%% remove any CARRIAGE RETURN.
            %%%%

            self.contents = strrep(fileread(file), char(13), '');
            self.contents = strrep(self.contents, newline, [' ', newline]);
            self.lineList = strsplit(self.contents, newline);
            self.lineListLen = length(self.lineList);
            iend = 1; % the line we are about to parse.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% determine the decoration symbol and read the banner.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%
            %%%% Here we assume the first non-empty line starts with the symbol character (likely ``%``).
            %%%%

            ibeg = self.skipNull(iend);
            failed = self.lineListLen < ibeg || length(self.lineList{ibeg}) < 4;
            if ~failed
                line = self.lineList{ibeg};
                failed = any(line(1) ~= line(2 : 4));
                if ~failed
                    iend = self.skipFull(ibeg + 1);
                    failed = self.lineListLen < iend;
                end
            end
            if ~failed
                self.dsym = line(1);
                self.banner = self.concat(self.lineList(ibeg : iend));
            else
                self.warn(line, "Failed to detect the decoration symbol used in the report file.");
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the ParaMonte MATLAB library interface specifications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.setup.library.interface = self.parseSection("interface.specifications");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the ParaMonte MATLAB library compiler version
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.setup.library.compiler.version = self.parseSection("compiler.version");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the ParaMonte MATLAB library compiler options
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.setup.library.compiler.options = self.parseSection("compiler.options");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the Runtime platform specifications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.setup.platform = self.parseSection("runtime.platform.specifications");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the simulation environment
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.setup.io = self.parseSection("simulation.environment");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the simulation specifications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %self.spec.base = self.parseSection("simulation.specifications.base");
            %if  self.method == "ParaDRAM" || self.method == "ParaDISE"
            %    self.spec.mcmc = self.parseSection("simulation.specifications.mcmc");
            %    self.spec.dram = self.parseSection("simulation.specifications.dram");
            %elseif self.method == "ParaNest"
            %    self.spec.nest = self.parseSection("simulation.specifications.nest");
            %end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% statistics: this must be always the last item to parse
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.checkpoint("parsing the simulation statistics...", 1);

            while true

                item = self.lineList{iend};
                if ~self.silent
                    self.spinner.spin(iend / self.lineListLen);
                end
                if  strcmp(item(1 : min(6, length(item))), 'stats.')

                    %%%%
                    %%%% First non-empty line.
                    %%%%

                    ibeg = self.skipNull(iend + 1);
                    failed = self.lineListLen < ibeg;
                    if ~failed
                        iend = self.skipFull(ibeg + 1);
                        failed = self.lineListLen < iend;
                        if ~failed
                            % Retrieve the value.
                            % The strategy is to read all values as table by dumping
                            % the value in a temporary file and reading it as a table.
                            tempfile = tempname();
                            fid = fopen(tempfile, "w");
                            fprintf(fid, self.concat(strtrim(self.lineList(ibeg : iend - 1))));
                            fclose(fid);
                            value = readtable(tempfile);
                            %self.(strtrim(item)).value = value;
                            eval(['self.', strtrim(item), '.value = value;']);
                        end
                    end

                    %%%%
                    %%%% Retrieve the description.
                    %%%%

                    ibeg = self.skipNull(iend + 1);
                    iend = self.skipFull(ibeg + 1);
                    if ~all(self.isdesc(self.lineList(ibeg : iend - 1)))
                        self.warn(string(ibeg) + "-" + string(iend), "A description record does not contain """ + string(self.prefix) + """");
                    else
                        % Concatenate lines and remove prefixes from each description lines.
                        desc = strrep(self.concat(strtrim(self.lineList(ibeg : iend - 1))), self.prefix, '');
                        eval(['self.', strtrim(item), '.description = desc;']);
                        %self.(strtrim(item)).description = desc;
                    end

                end % new item found

                iend = iend + 1;
                if  self.lineListLen < iend
                    break;
                end

            end % while

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.checkpoint();

            %%%%
            %%%% Add the parallel scaling visualizations.
            %%%%

            silent_kws = {"silent", self.silent};

            try
                %self.vis = struct();
                %self.vis.parallelism = struct();
                %self.vis.parallelism = struct();
                for parcond = ["optimal", "perfect"]
                    %self.vis.parallelism.(parcond) = struct();
                    %self.vis.parallelism.(parcond).scaling = struct();
                    %self.vis.parallelism.(parcond).scaling.strong = struct();
                    %self.vis.parallelism.(parcond).scaling.strong.lineScatter = ...
                    self.stats.parallelism.(parcond).scaling.strong.vis = struct();
                    self.stats.parallelism.(parcond).scaling.strong.vis.desc = ...
                    "This component is added in the post-processing phase to facilitate quick " + ...
                    "visualization of the parallel scaling behavior of the sampler under the " + parcond + " parallel conditions.";
                    self.stats.parallelism.(parcond).scaling.strong.vis.lineScatter = ...
                    pm.vis.PlotLineScatter  ( self.stats.parallelism.(parcond).scaling.strong.value ...
                                            , "colx", 1, "coly", 2, "colc", 2, "axes", {"xscale", "log"} ...
                                            , "scatter", {"size", 10, "color", pm.vis.color.rgb("red")} ...
                                            ..., "colormap", {"enabled", false, "map", "autumn"} ...
                                            , "plot", {"linewidth", 2.5} ...
                                            , silent_kws{:} ...
                                            );
                end
            catch me
                warning ( newline ...
                        + "Failed to create the visualizations for the parallelization scaling behavior of the sampler." + newline ...
                        + "It is possible that the component ``stats.parallelism.(parcond).scaling.strong.value`` of the " + newline ...
                        + "output ``report`` object of class ``pm.sampling.FileContentsReport`` does not exist, where " + newline ...
                        + "the string `parcond` can be either ``optimal`` or ``perfect``." + newline ...
                        + "Here is the error message:" + newline ...
                        + newline ...
                        + string(me.identifier) + newline + string(me.message) + newline ...
                        + newline ...
                        );
            end

            %%%%
            %%%% Add the parallel process contribution visualizations.
            %%%%

            try
                %%%%
                %%%% Set up the Cyclic Geometric fit to parallel processes contributions.
                %%%%

                cgfit = @(iproc, successProb, normFac, processCount) successProb .* normFac .* (1 - successProb).^(iproc - 1) ./ (1 - (1 - successProb).^processCount);

                %%%%
                %%%% Set up the data frame containing the processes contributions and its Cyclic Geometric fit to parallel processes contributions.
                %%%%

                dfpc = self.stats.parallelism.current.process.contribution.count.value;
                dfpc.cyclicGeometricFit = cgfit ( dfpc.processID ...
                                                , self.stats.parallelism.current.process.contribution.fit.value.successProb ...
                                                , self.stats.parallelism.current.process.contribution.fit.value.normFac ...
                                                , self.stats.parallelism.current.process.count.value.Var1 ...
                                                );
                self.stats.parallelism.current.process.contribution.vis = struct();
                self.stats.parallelism.current.process.contribution.vis.desc = ...
                "This component is added in the post-processing phase to facilitate quick " + ...
                "visualization of the contributions of the parallel processes to the final chain.";
                self.stats.parallelism.current.process.contribution.vis.line = ...
                pm.vis.PlotLine ( dfpc ...
                                , "colx", "processID", "coly", 2:3, "axes", {"xscale", "log", "yscale", "log"} ...
                                , "ylabel", {"txt", "Accepted-Sample Contribution"} ...
                                , "xlabel", {"txt", "Process ID"} ...
                                , "colormap", {"enabled", false} ...
                                , "legend", {"enabled", true} ...
                                , "plot", {"linewidth", 2.5} ...
                                , silent_kws{:} ...
                                );
            catch me
                warning ( newline ...
                        + "Failed to create the visualizations for the parallel processes contributions and their Cyclic Geomtric fit." + newline ...
                        + "It is possible that the component ``stats.parallelism.current.process.contribution`` of the " + newline ...
                        + "output ``report`` object of class ``pm.sampling.FileContentsReport`` does not exist, or " + newline ...
                        + "its sub-components are missing." + newline ...
                        + "Here is the error message:" + newline ...
                        + newline ...
                        + string(me.identifier) + newline + string(me.message) + newline ...
                        + newline ...
                        );
            end

        end % constructor

    end % methods(Access = public)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Hidden, Static)

        %function reportParseFailure(topic)
        %    topic = string(strtrim(strrep(strrep(topic, newline, ' '), char(13), ' ')));
        %    warning ( newline ...
        %            + "Failed to parse the record """ + topic + """. " + newline ...
        %            + "The structure of the report file appears to have been compromised. Skipping..." + newline ...
        %            + newline ...
        %            );
        %end
        %
        %function reportMissingValue(topic)
        %    topic = string(strtrim(strrep(strrep(topic, newline, ' '), char(13), ' ')));
        %    warning ( newline ...
        %            + "Failed to parse the value corresponding to the record """ + topic + """. " + newline ...
        %            + "The structure of the report file appears to have been compromised. Skipping..." + newline ...
        %            + newline ...
        %            );
        %end

        function str = concat(lines)
            % Concatenate the elements of an input cell array by
            % the new line character and return the result as a single string.
            str = string(join(lines, newline));
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Hidden)

        %function result = isSectionHeader(self, record)
        %    result = self.dsymLen < length(record) && contains(record(1 : self.dsymLen), self.dsym);
        %end
        %
        %function stopBeforeNextSectionHeader(self)
        %    % NOTE: if it is already within a section header, it will skip it and find the next section header.
        %    record = self.lineList{iend};
        %    isCurrentSectionHeader = self.dsymLen < length(record) && contains(record(1 : self.dsymLen), self.dsym);
        %    while true
        %        if iend == self.lineListLen
        %            break;
        %        end
        %        record = self.lineList{iend};
        %        if  self.dsymLen < length(record) && contains(record(1 : self.dsymLen), self.dsym)
        %            if isCurrentSectionHeader
        %                iend = iend + 1;
        %                continue;
        %            else
        %                iend = iend - 1;
        %                break;
        %            end
        %        else
        %            isCurrentSectionHeader = false;
        %            iend = iend + 1;
        %            continue;
        %        end
        %    end
        %end
        %
        %function skipCurrentSectionHeader(self)
        %    while true
        %        record = self.lineList{iend};
        %        if  self.dsymLen < length(record) && contains(record(1 : self.dsymLen), self.dsym)
        %            iend = iend + 1;
        %            if iend == self.lineListLen
        %                break;
        %            end
        %            continue;
        %        else
        %            break;
        %        end
        %    end
        %end
        %
        %function section = parseSection(self, topic)
        %    ilineLastSuccess = iend;
        %    topic = lower(topic);
        %    section = [];
        %    topicFound = false;
        %    while true
        %        iend = iend + 1;
        %        if  self.lineListLen <= iend
        %            break;
        %        end
        %        record = lower(self.lineList{iend});
        %        if contains(record, topic)
        %            topicFound = true;
        %            break;
        %        else
        %            continue;
        %        end
        %    end
        %    if  topicFound
        %        self.skipCurrentSectionHeader();
        %        lineStart = iend;
        %        self.stopBeforeNextSectionHeader();
        %        section = self.concat(lineStart, iend);
        %    else
        %        self.reportParseFailure(topic);
        %        iend = ilineLastSuccess;
        %    end
        %end

        function result = isdesc(self, record)
            result = contains(record, self.prefix);
        end

        function iend = skipNull(self, ibeg)
            % skip the empty lines after the current
            % until encountering the first non-empty line.
            iend = ibeg;
            while true
                if  self.lineListLen < iend
                    break;
                elseif ~isempty(strtrim(self.lineList{iend}))
                    break;
                end
                iend = iend + 1;
            end
        end

        function iend = skipFull(self, ibeg)
            % skip the nonempty lines after the current
            % until encountering the first empty line.
            iend = ibeg;
            while true
                if  self.lineListLen < iend
                    break;
                elseif isempty(strtrim(self.lineList{iend}))
                    break;
                end
                iend = iend + 1;
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end