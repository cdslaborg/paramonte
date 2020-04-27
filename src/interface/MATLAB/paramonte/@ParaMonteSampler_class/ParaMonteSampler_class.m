classdef ParaMonteSampler_class < dynamicprops
%   Base class for the ParaMonte sampler routines

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        mpiEnabled = false;
        inputFile = [];
        buildMode = "release";
        spec = [];
    end

    properties (Access = public, Hidden)
        Err = Err_class();
        methodName = "";
        objectName = [];
        ndim = [];
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)
        [chainList] = readChain(self,varargin)
        [sampleList] = readSample(self,varargin)
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)
        % These methods have been implemented to override the default 'handle' class methods, 
        % so that they won't pop-up after pressing 'Tab' button.
        function addlistener    (); end
        function delete         (); end
        function findobj        (); end
        function findprop       (); end
        function valid          (); end
        function listener       (); end
        function notify         (); end
        function name = getMyName(self); name = inputname(1); end
        fileList = getFileList(self,file,fileType)
        result = genOutputFileName(self)
        namelist = getInputFile(self)
        outputList = readOutput(self,file,delimiter,fileType)
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef ParaMonteSamplerRoutine
