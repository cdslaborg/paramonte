function out = parseArgs(self,varargin)

    vararginLen = length(varargin);
    if isa(self,"struct")
        selfProperties = fieldnames(self);
        selfPropertiesLen = length(selfProperties);
    else
        selfProperties = properties(self);
        selfPropertiesLen = length(selfProperties);
    end

    for i = 1:2:vararginLen
        propertyDoesNotExist = true;
        vararginString = string(varargin{i});
        for ip = 1:selfPropertiesLen
            if strcmp(vararginString,string(selfProperties(ip)))
                propertyDoesNotExist = false;
                if i < vararginLen
                    if isa(self.(selfProperties{ip}),"struct") && isa(varargin{i+1},"cell")
                        self.(selfProperties{ip}) = parseArgs( self.(selfProperties{ip}) , varargin{i+1}{:} );
                    else
                        self.(selfProperties{ip}) = varargin{i+1};
                    end
                else
                    error("The corresponding value for the property """ + string(selfProperties{ip}) + """ is missing as input argument.");
                end
                break;
            end
        end
        if propertyDoesNotExist
            error("The requested property """ + string(varargin{i}) + """ does not exist.");
        end
    end

    out = self;

end
