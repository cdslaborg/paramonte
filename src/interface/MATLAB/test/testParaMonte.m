clc;
clear all;
close all;
clear classes;
format compact; format long;
pmlibRootDir = '../'; % set this path to the paramonte library root dir
addpath( genpath(pmlibRootDir) ) 
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath); % Change working directory to source code directory.
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

%-----------------------------------------------------------------------------------------------------------------------------------------------------------
%pm = paramonte("matlab");
pm = paramonte(); % "matlab"
pmpd = pm.ParaDRAM();
%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%pmpd.inputFile = './paramonte.in';
pmpd.spec.chainSize = 10000;
pmpd.spec.adaptiveUpdateCount = 0;

%pmpd.runSampler(2,@getLogFunc);
%pmpd.readChain("D:\Dropbox\Projects\20180101_ParaMonte\git\src\interface\MATLAB\test\ParaDRAM_run_180420_213248_768");

%pmpd.spec.chainSize                             = 10000;
%pmpd.spec.outputFileName                        = ["./out/temp/"];                % Works
%pmpd.runSampler(4,@getLogFunc);
%pmpd.readChain("ParaDRAM_run_")

pmpd.spec.chainSize                             = 10000;
pmpd.spec.adaptiveUpdateCount                   = 0;
%pmpd.runSampler(20,@getLogFunc)
pmpd.spec.outputFileName = "ParaDRAM_run_180420_213248_768";
pmpd.readMarkovChain();
pmpd.readSample();
pmpd.readChain(); %"ParaDRAM_run_180420_151344_607";%"D:\Dropbox\Projects\20180101_ParaMonte\git\src\interface\MATLAB\test\ParaDRAM_run_160420_023054_530");
return
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaMonte variables...
% pmpd.spec.sampleSize                            = 1000;                         % Works
% pmpd.spec.randomSeed                            = 7;                            % Works
pmpd.spec.description                           = "Hi there";                   % Works
pmpd.spec.outputFileName                        = "./out/temp/";                % Works
pmpd.spec.outputDelimiter                       = "|";                          % Works
pmpd.spec.chainFileFormat                       = 'verbose';                    % Works
pmpd.spec.variableNameList                      = ["Variable-X", "Variable-Y"]; % Works
% pmpd.spec.restartFileFormat                     = [];                           % Not implemented properly yet.
pmpd.spec.outputColumnWidth                     = 13;                           % Works
pmpd.spec.outputRealPrecision                   = 6;                            % Works
pmpd.spec.silentModeRequested                   = 0;                            % Works
pmpd.spec.domainLowerLimitVec                   = [-4,-4];                      % Works
pmpd.spec.domainUpperLimitVec                   = [+4,+4];                      % Works
pmpd.spec.progressReportPeriod                  = 100;                          % Works
pmpd.spec.targetAcceptanceRate                  = 0.5;                          % Works
% pmpd.spec.maxNumDomainCheckToWarn               = [];
% pmpd.spec.maxNumDomainCheckToStop               = [];
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaMCMC variables...
pmpd.spec.chainSize                             = 5000;                         % Works
pmpd.spec.startPointVec                         = [0.6,1.2];                    % Works
pmpd.spec.sampleRefinementCount                 = 1;                            % Works
% pmpd.spec.sampleRefinementMethod                = "someRandomName";             % Works
pmpd.spec.randomStartPointRequested             = 1;                            % Works
pmpd.spec.randomStartPointDomainLowerLimitVec   = [0.5, 1.0];                   % Works
pmpd.spec.randomStartPointDomainUpperLimitVec   = [1.0, 1.5];                   % Works
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaDRAM variables...
pmpd.spec.scaleFactor                           = "1.2 * gelman";               % Works
pmpd.spec.proposalModel                         = "normal";                     % Works
pmpd.spec.proposalStartCovMat                   = [0.5, 0.2; 0.1, 0.3];         % Works
pmpd.spec.proposalStartCorMat                   = [0.4, 0.1; 0.2, 0.3];         % Works
pmpd.spec.proposalStartStdVec                   = [1, 1];                       % Works
pmpd.spec.adaptiveUpdateCount                   = 2;                            % Works
pmpd.spec.adaptiveUpdatePeriod                  = 25;                           % Works
pmpd.spec.greedyAdaptationCount                 = 2;                            % Works
pmpd.spec.delayedRejectionCount                 = 2;                            % Works
pmpd.spec.burninAdaptationMeasure               = 0.5;                          % Works
pmpd.spec.delayedRejectionScaleFactorVec        = [3, 4];                       % Works
%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndim = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmpd.runSampler(ndim, @getLogFunc);
fclose('all');

for i = 1 : ndim

    figure;
    histfit(pmpd.Chain.State(i,:));
    
    figure;
    plot(pmpd.Chain.State(i,:));
    set(gca,'xscale','log')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system(pmpd.LogFile.Path.modified);
system(pmpd.TimeFile.Path.modified);
system(pmpd.ChainFile.Path.modified);
system(pmpd.SampleFile.Path.modified);
