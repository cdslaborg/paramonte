function parseArgs(self,varargin)

    vararginLen = length(varargin);
    for i = 1:2:vararginLen
        propertyDoesNotExist = true;
        selfProperties = properties(self);
        selfPropertiesLen = length(selfProperties);
        for ip = 1:selfPropertiesLen
            if strcmp(string(varargin{i}),string(selfProperties{ip}))
                propertyDoesNotExist = false;
                if i < vararginLen
                    self.(selfProperties{ip}) = varargin{i+1};
                else
                    error("The corresponding value for the property """ + string(selfProperties{ip}) + """ is missing as input argument.");
                end
                break;
            end
        end
        if propertyDoesNotExist
            error("The requested the property """ + string(varargin{i}) + """ does not exist.");
        end
    end

end
