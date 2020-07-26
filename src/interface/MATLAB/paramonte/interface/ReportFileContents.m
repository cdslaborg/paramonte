%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it 
%   under the terms of the GNU Lesser General Public License as published 
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see, 
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms, 
%   if you use any parts of this library for any purposes, 
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This is the ReportFileContents class for generating instances 
%   of ParaMonte output report file contents. The ParaMonte readReport method
%   returns an object or a list of objects of class ReportFileContents.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef ReportFileContents < OutputFileContents

    properties(Access = public)
        contents = [];
        setup = struct();
        stats = struct();
        spec = struct();
        lineList = [];
    end

    properties(Hidden)
        lineListLen = [];
        indentLen = 8; % indent length of the records
        dsym = '****'; % decoration symbol
        lineCounter;
        dsymLen;
        prefix;
        newline = char(10);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = ReportFileContents  ( file ...
                                            , methodName ...
                                            , mpiEnabled ...
                                            , Err ...
                                            )

            self = self@OutputFileContents(file,methodName,mpiEnabled,Err);
            self.timer.tic();

            self.dsymLen = length(self.dsym); % decoration symbol
            self.prefix = convertStringsToChars(self.methodName + " - NOTE:");

            self.contents = strrep(fileread(file),char(13),'');
            %if ispc
            %    self.contents = strrep(fileread(file),char(13),'');
            %else
            %    self.contents = fileread(file);
            %end
            self.contents = strrep(self.contents,self.newline,[' ',self.newline]);
            self.lineList = strsplit(self.contents,self.newline); % strtrim()
            self.lineListLen = length(self.lineList);

            self.updateUser("parsing the report file contents...");

            self.lineCounter  = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the banner
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lineStartFound = false;
            while true
                self.lineCounter = self.lineCounter + 1; if self.lineCounter>=self.lineListLen; break; end
                record = self.lineList{self.lineCounter};
                if lineStartFound
                    if ~contains(record,self.dsym)
                        self.lineCounter = self.lineCounter - 1;
                        break
                    end
                else
                    if contains(record,self.dsym)
                        lineStartFound = true;
                        lineStart = self.lineCounter;
                    end
                end
            end
            if lineStartFound
                self.setup.library.banner = self.concat(lineStart,self.lineCounter);
            else
                self.reportParseFailure("ParaMonte banner");
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the ParaMonte library interface specifications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.setup.library.interface = self.parseSection("ParaMonte library interface specifications");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the Runtime platform specifications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.setup.platform = self.parseSection("Runtime platform specifications");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the simulation environment
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.setup.io = self.parseSection("simulation environment");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% read the simulation environment
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.spec = self.parseSection("simulation specifications");

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% statistics: this must be always the last item to parse
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.parseStats();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.updateUser([]);

        end % constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %function helpme(self,varargin)
        %    %
        %    %   Open the documentation for the input object's name in string format, otherwise, 
        %    %   open the documentation page for the class of the object owning the helpme() method.
        %    %
        %    %   Parameters
        %    %   ----------
        %    %
        %    %       This function takes at most one string argument, 
        %    %       which is the name of the object for which help is needed.
        %    %
        %    %   Returns
        %    %   -------
        %    %
        %    %       None. 
        %    %
        %    %   Example
        %    %   -------
        %    %
        %    %       helpme("plot")
        %    %
        %    methodNotFound = true;
        %    if nargin==2
        %        if strcmpi(varargin{1},"reset")
        %            cmd = "doc self.resetPlot";
        %            methodNotFound = false;
        %        else
        %            methodList = ["plot","helpme"];
        %            for method = methodList
        %                if strcmpi(varargin{1},method)
        %                    methodNotFound = false;
        %                    cmd = "doc self." + method;
        %                end
        %            end
        %        end
        %    elseif nargin~=1
        %        error("The helpme() method takes at most one argument that must be string.");
        %    end
        %    if methodNotFound
        %        cmd = "doc self";
        %    end
        %    eval(cmd);
        %end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = public)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reportParseFailure(self,topic)
            topic = string(strtrim(strrep(strrep(topic,self.newline,' '),char(13),' ')));
            self.Err.msg    = "Failed to parse the record """ + topic + """. " ...
                            + "The structure of the report file appears to have been compromised. Skipping... "; 
            self.Err.warn();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function reportMissingValue(self,topic)
            topic = string(strtrim(strrep(strrep(topic,self.newline,' '),char(13),' ')));
            self.Err.msg    = "Failed to parse the value corresponding to the record """ + topic + """. " ...
                            + "The structure of the report file appears to have been compromised. Skipping... "; 
            self.Err.warn();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function result = isSectionHeader(self,record)
            if length(record)>self.dsymLen && contains(record(1:self.dsymLen),self.dsym)
                result = true;
            else
                result = false;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function stopBeforeNextSectionHeader(self)
            % NOTE: if it is already within a section header, it will skip it and find the next section header.
            record = self.lineList{self.lineCounter};
            isCurrentSectionHeader = length(record)>self.dsymLen && contains(record(1:self.dsymLen),self.dsym);
            while true
                if self.lineCounter==self.lineListLen; break; end
                record = self.lineList{self.lineCounter};
                if length(record)>self.dsymLen && contains(record(1:self.dsymLen),self.dsym)
                    if isCurrentSectionHeader
                        self.lineCounter = self.lineCounter + 1;
                        continue;
                    end
                    self.lineCounter = self.lineCounter - 1;
                    break;
                else
                    isCurrentSectionHeader = false;
                    self.lineCounter = self.lineCounter + 1;
                    continue;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function skipCurrentSectionHeader(self)
            while true
                record = self.lineList{self.lineCounter};
                if length(record)>self.dsymLen && contains(record(1:self.dsymLen),self.dsym)
                    self.lineCounter = self.lineCounter + 1; if self.lineCounter==self.lineListLen; break; end
                    continue;
                else
                    break;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function section = parseSection(self,topic)
            lineCounterLastSuccess = self.lineCounter;
            section = [];
            topicFound = false;
            while true
                self.lineCounter = self.lineCounter + 1;
                if self.lineCounter>=self.lineListLen; break; end
                record = self.lineList{self.lineCounter};
                if contains(record,topic)
                    topicFound = true;
                    break;
                else
                    continue;
                end
            end
            if topicFound
                self.skipCurrentSectionHeader();
                lineStart = self.lineCounter;
                self.stopBeforeNextSectionHeader();
                section = self.concat(lineStart,self.lineCounter);
            else
                self.reportParseFailure(topic);
                self.lineCounter = lineCounterLastSuccess;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function concatatedString = concat(self,lineStart,lineCounter)
            concatatedString = self.lineList(lineStart:lineCounter);
            %concatatedString = concatatedString(~cellfun('isempty',concatatedString));
            concatatedString = string(join(concatatedString,self.newline));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function parseStats(self)

            self.stats = struct();
            if strcmp(self.methodName,"ParaDRAM") || strcmp(self.methodName,"MatDRAM")

                while true

                    self.lineCounter = self.lineCounter + 1; if self.lineCounter>self.lineListLen; break; end
                    item = self.lineList{self.lineCounter};

                    if self.isstats(item)

                        % parse the value

                        valueIsNumeric = true;
                        valueFound = false;
                        descFound = false;
                        value = '';
                        desc = '';
                        while true

                            self.lineCounter = self.lineCounter + 1; if self.lineCounter>self.lineListLen; break; end
                            record = strrep(strrep(self.lineList{self.lineCounter},self.newline,' '), char(13), ' '); % replace newline with space;

                            % check the record is not another item or desc is not found before value.

                            if self.isstats(record) %&& ~(descFound && valueFound)
                                %self.reportParseFailure(item);
                                self.lineCounter = self.lineCounter - 1;
                                break
                            end

                            % parse the value/description

                            if length(record)>self.indentLen
                                recordIsDesc = self.isdesc(record);
                                if recordIsDesc
                                    if ~valueFound
                                        self.reportMissingValue(item);
                                        %self.lineCounter = self.lineCounter - 1;
                                        break
                                    end
                                    descFound = true;
                                    desc = [ desc, ' ', strtrim( strrep(self.lineList{self.lineCounter},self.prefix,' ') ) ]; % remove prefix, trim, append.
                                elseif ~descFound
                                    valueFound = true;
                                    value = [value, record];
                                    if isempty(str2num(record))
                                        valueIsNumeric = false;
                                    end
                                end
                            end

                        end

                        if valueFound && descFound
                            if valueIsNumeric
                                value = ['[',value,']'];
                                %value = strsplit(join(strtrim(self.lineList(valueStart:valueEnd)),' '));
                                %value = value(~cellfun('isempty',value));
                                eval(['self.',strtrim(item),'.value=',value,';']);
                                eval(['self.',strtrim(item),'.description="',strtrim(desc),'";']);
                                %valueStart = self.lineCounter;
                            %else
                            %    value = [];
                            %    for i = valueStart:valueEnd
                            %        value = [value, self.lineList{i}, self.newline];
                            %    end
                            %    value = ['''', value, ''''];
                            %    disp(['self.',strtrim(item),'.value=''',value,''';']);
                            %    eval(['self.',strtrim(item),'.value=''',value,''';']);
                            end
                        end

                    end % new item found

                end % while

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function result = isstats(self,record)
            result = length(record)>5 && strcmp(record(1:6),'stats.');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function result = isdesc(self,record)
            result = contains(record,self.prefix);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Access = public)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef ReportFileContents < handle