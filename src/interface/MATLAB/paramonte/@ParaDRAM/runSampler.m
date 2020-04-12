function runSampler(ndim,getLogFunc)
%   Run ParaDRAM sampler and return nothing.
%   
%   Parameters
%   ----------
%       ndim
%           integer representing the number of dimensions of the
%           domain of the user's objective function getLogFunc().
%           It must be a positive integer.
%       getLogFunc()
%           represents the user's objective function to be sampled,
%           which must take a single input argument of type numpy
%           float64 array of length ndim and must return the
%           natural logarithm of the objective function.
%   
%   Returns
%   -------
%       None
end