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
