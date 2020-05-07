function runSampler(self,ndim,getLogFunc,varargin)
    self.objectName = inputname(1);
    runSampler@ParaMonteSampler(self,ndim,getLogFunc,varargin{:})
end