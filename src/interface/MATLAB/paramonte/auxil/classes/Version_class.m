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
classdef Version_class

    properties(Access=public)
    end

    properties(Access=public, Hidden)
        versionList = ["interface","kernel"];
        versionPath
        versionType
        versionSave
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Version_class(versionPath,versionType)
            self.versionPath = string(versionPath);
            self.versionType = string(versionType);
            self.versionSave = [];
            self.checkVersionType();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function version = get(self)
            version = "ParaMonte MATLAB " + string([upper(self.versionType{1}(1)),self.versionType{1}(2:end)]) + " Version " + self.dump();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function numericVersion = dump(self)
            for i = 1:length(self.versionList)
                versionType = self.versionList(i);
                if strcmpi(self.versionType,versionType)
                    if isempty(self.versionSave)
                        versionFileName = ".VERSION_" + upper(self.versionList(i));
                        break;
                    else
                        numericVersion = self.versionSave;
                        return;
                    end
                end
            end
            versionFilePath = convertStringsToChars(fullfile(self.versionPath, versionFileName));
            try
                fid = fopen(versionFilePath);
                numericVersion = string(fgetl(fid));
                fclose(fid);
            catch
                numericVersion = "UNKNOWN";
            end
            self.versionSave = numericVersion;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function thisVersionType = checkVersionType(self)

            versionTypeNotFound = true;
            for versionType = self.versionList
                if strcmpi(versionType,self.versionType)
                    versionTypeNotFound = false;
                    break;
                end
            end
            if versionTypeNotFound
                error("The input versionType is not a valid recognized version type. Possible values: " + join(versionList," "));
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
