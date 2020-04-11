clc;
clear all;
close all;
clear classes;
format compact; format long;
addpath(genpath('..\.')) % lib codes
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath); % Change working directory to source code directory.
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

%-----------------------------------------------------------------------------------------------------------------------------------------------------------
                pd = ParaDRAM();
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaMonte variables...
% pd.spec.sampleSize                              = 1000;                         % Works
pd.spec.randomSeed                              = 7;                            % Works
pd.spec.description                             = "Hi there";                   % Works
pd.spec.outputFileName                          = "./out/temp/";                % Works
pd.spec.outputDelimiter                         = "|";                          % Works
pd.spec.chainFileFormat                         = 'verbose';                    % Works
pd.spec.variableNameList                        = ["Variable-X", "Variable-Y"]; % Works
% pd.spec.restartFileFormat                       = [];                         % Not implemented properly yet.
pd.spec.outputColumnWidth                       = 13;                           % Works
pd.spec.outputRealPrecision                     = 6;                            % Works
pd.spec.silentModeRequested                     = 0;                            % Works
pd.spec.domainLowerLimitVec                     = [-4,-4];                      % Works
pd.spec.domainUpperLimitVec                     = [+4,+4];                      % Works
pd.spec.progressReportPeriod                    = 100;                          % Works
pd.spec.targetAcceptanceRate                    = 0.5;                          % Works
% pd.spec.maxNumDomainCheckToWarn                 = [];
% pd.spec.maxNumDomainCheckToStop                 = [];
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaMCMC variables...
pd.spec.chainSize                               = 5000;                         % Works
pd.spec.startPointVec                           = [0.6,1.2];                    % Works
pd.spec.sampleRefinementCount                   = 1;                            % Works
% pd.spec.sampleRefinementMethod                  = "SomeRandomName";             % Works
pd.spec.randomStartPointRequested               = 1;                            % Works
pd.spec.randomStartPointDomainLowerLimitVec     = [0.5, 1.0];                   % Works
pd.spec.randomStartPointDomainUpperLimitVec     = [1.0, 1.5];                   % Works
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaDRAM variables...
% pd.spec.scaleFactor                             = "2.5";                        % Problem: Doesn't work for '1', "1", '2.5', "2.5", '2.5*Gelman' or "2.5*Gelman"
pd.spec.proposalModel                           = "uniform";                    % Works
pd.spec.proposalStartCovMat                     = [0.5, 0.2; 0.1, 0.3];         % Works
pd.spec.proposalStartCorMat                     = [0.4, 0.1; 0.2, 0.3];         % Works
pd.spec.proposalStartStdVec                     = [1, 1];                       % Works
pd.spec.adaptiveUpdateCount                     = 2;                            % Works
pd.spec.adaptiveUpdatePeriod                    = 25;                           % Works
pd.spec.greedyAdaptationCount                   = 2;                            % Works
pd.spec.delayedRejectionCount                   = 2;                            % Works
pd.spec.burninAdaptationMeasure                 = 0.5;                          % Works
pd.spec.delayedRejectionScaleFactorVec          = [3, 4];                       % Works
%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndim = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pd.runSampler(ndim, @getLogFunc);
fclose('all');

for i = 1 : ndim

    figure;
    histfit(pd.Chain.State(i,:));
    
    figure;
    plot(pd.Chain.State(i,:));
    set(gca,'xscale','log')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(pd.LogFile.Path.modified);
system(pd.TimeFile.Path.modified);
system(pd.ChainFile.Path.modified);
system(pd.SampleFile.Path.modified);
