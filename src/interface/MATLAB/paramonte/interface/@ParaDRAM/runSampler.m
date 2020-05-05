function runSampler(self,ndim,getLogFunc,varargin)
    %functions(getLogFunc)
    %global getLogFuncHandle
    %getLogFuncHandle = getLogFuncHandleArg;
    %return
    runSampler@ParaMonteSampler(self,ndim,getLogFunc,varargin{:})
end