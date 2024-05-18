clc;
clear all;
close all;
clear classes;
format compact; format long;
pmdir = './'; addpath(pmdir); % change this path to the ParaMonte library root directory
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

sampler = pm.sampling.Paradram();
sampler.spec.randomSeed = 284651; % make sampling reproducible.
%sampler.spec.outputChainSize = 100000; % Use a small chain size for illustration.
sampler.spec.parallelismNumThread = []; % Use these many parallel threads.
sampler.run(@getLogFunc, 2);

function logFunc = getLogFunc(state)
    %
    % Return the negative natural logarithm of the 2-dimensional Himmelblau function.
    % Reference: https://en.wikipedia.org/wiki/Himmelblau%27s_function
    %
    logFunc = -log((state(1)^2 + state(2) - 11)^2 + (state(1) + state(2)^2 - 7)^2 + 0.1);
    if false
        nsim = 100000;
        for i = 1 : nsim
            logFunc = logFunc - log((state(1)^2 + state(2) - 11)^2 + (state(1) + state(2)^2 - 7)^2 + 0.1);
        end
        logFunc = logFunc / (nsim + 1);
    end
end
