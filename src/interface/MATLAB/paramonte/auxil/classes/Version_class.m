classdef Version_class

    properties(Access=public)
        interface
        kernel
    end

    properties(Access=public, Hidden)
        root
        savedStruct
        versionList = ["interface","kernel"];
        versionFile = [".VERSION",".VERSION_KERNEL"];
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = Version_class(rootPath)
            self.root = rootPath;
            self.savedStruct = struct();
            for version = self.versionList
                self.savedStruct.(version) = [];
                self.(version) = self.get(version);
            end
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function version = get(self,varargin)
            this = self.checkVersionType("get",varargin{:});
            version = "ParaMonte MATLAB " + string([upper(this{1}(1)),this{1}(2:end)]) + " Version " + self.dump(this);
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function numericVersion = dump(self,varargin)
            this = self.checkVersionType("dump",varargin{:});
            for i = 1:length(self.versionList)
                version = self.versionList(i);
                if strcmpi(this,version)
                    if isempty(self.savedStruct.(this))
                        versionFileName = self.versionFile(i);
                        break;
                    else
                        numericVersion = self.savedStruct.(this);
                        return;
                    end
                end
            end
            versionFilePath = convertStringsToChars(fullfile(self.root, versionFileName));
            try
                versionFileID = fopen(versionFilePath);
                numericVersion = fgetl(versionFileID);
            catch
                numericVersion = "UNKNOWN";
            end
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public, Hidden)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function thisVersionType = checkVersionType(self,versionClassMethodName,varargin)

            invalidInputProvided = false;
            if nargin==2
                thisVersionType = self.versionList(1);
            elseif nargin==3
                thisVersionType = varargin{1};
            else
                invalidInputProvided = true;
            end

            if ~invalidInputProvided
                if getVecLen(thisVersionType)
                    invalidString = ~(isa(thisVersionType,"string") && length(thisVersionType)==1 && any(strcmpi(thisVersionType,self.versionList)));
                    invalidChar = ~(isa(thisVersionType,"char") && any(strcmpi(string(thisVersionType),self.versionList)));
                    if invalidString && invalidChar; invalidInputProvided = true; end
                else
                    thisVersionType = self.versionList(1);
                end
            end

            if invalidInputProvided
                error   ( newline ...
                        + "The " + versionClassMethodName + "() method of Version_class() only takes one input string argument with two possible values:" + newline + newline ...
                        + "    " + versionClassMethodName + "(""" + versionList(1) + """)" + newline ...
                        + "    " + versionClassMethodName + "(""" + versionList(2) + """)" + newline ...
                        );
            else
                thisVersionType = string(thisVersionType);
            end

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
