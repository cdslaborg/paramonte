%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                                                                            %%%%
%%%%    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           %%%%
%%%%                                                                                                                            %%%%
%%%%    Copyright (C) 2012-present, The Computational Data Science Lab                                                          %%%%
%%%%                                                                                                                            %%%%
%%%%    This file is part of the ParaMonte library.                                                                             %%%%
%%%%                                                                                                                            %%%%
%%%%    LICENSE                                                                                                                 %%%%
%%%%                                                                                                                            %%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          %%%%
%%%%                                                                                                                            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear MATLAB space

clc;
clear all;
close all;
clear classes;
format compact; format long;

% set path to the ParaMonte MATLAB library

pmlibRootDir = './'; % set this path to the paramonte library root dir
addpath(genpath(pmlibRootDir));

% change MATLAB's working directory to the folder containing this script

filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir); % Change working directory to source code directory.
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

% define the objective function

NDIM = 2; % number of dimensions of the domain of the distribution

% create a ParaMonte object

pm = paramonte(); % use paramonte("matlab") to use the MATLAB kernel instead of interface

% create a ParaDRAM simulation object

pmpd = pm.ParaDRAM();

pmpd.mpiEnabled = false;
pmpd.bldtype = "release";
pmpd.spec.randomSeed = 31731;
pmpd.spec.outputFileName = "./RosenBrock/";
pmpd.spec.sampleRefinementCount = 1;
pmpd.spec.outputSeparator = ",";
pmpd.spec.chainSize = 15000;
pmpd.spec.sampleSize = -100;
pmpd.spec.targetAcceptanceRate = [0.01,0.3];
pmpd.spec.proposalStd = [2,2];

getNegLogRosen = @(point) -log( 100.0 * (point(2)-point(1)^2)^2 + (point(1)-1.0)^2 );

pmpd.runSampler ( NDIM              ... number of dimensions of the objective function
                , getNegLogRosen    ... the objective function
                );

function logFunc = getLogFunc(point)
    logFunc = -log( 100.0 * (point(2)-point(1)^2)^2 + (point(1)-1.0)^2 );
end
