function run(self, getLogFunc, ndim)
    %
    %   Run the sampler and return nothing.
    %
    %   Parameters
    %   ----------
    %
    %       getLogFunc()
    %
    %           The input MATLAB function handle or anonymous (lambda) function
    %           containing the implementation of the objective function to be sampled.
    %           This user-specified function must have the following interface,
    %
    %               function logFunc = getLogFunc(state)
    %               end
    %
    %           where,
    %
    %               1.  the input argument ``state`` is a vector of type MATLAB ``double``
    %                   of size ``ndim`` representing a single point from within the ``ndim``
    %                   dimensional domain of the mathematical object function to be explored.
    %
    %               2.  the output argument `logFunc` is a scalar of the same type as the
    %                   input ``state`` containing the natural logarithm of the objective
    %                   function at the specified input ``state`` within its domain.
    %
    %       ndim
    %
    %           The input scalar positive-valued whole-number representing the number of dimensions
    %           of the domain of the user-specified objective function in the input ``getLogFunc()``.
    %
    %   Returns
    %   -------
    %
    %       None
    %
    %   Interface
    %   ---------
    %
    %       sampler = pm.sampling.Paradram();
    %       sampler.run(getLogFunc, ndim);
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if nargin < 3
        help("pm.sampling.Paradram.run");
        error   ( newline ...
                + "The `run()` method of the " + self.method + " sampler" + newline ...
                + "requires at least two input arguments: ``getLogFunc()`` and ``ndim``" + newline ...
                + "For more information on the input arguments," + newline ...
                + "see the method documentation shown above." + newline ...
                + newline ...
                );
    end
    self.name = inputname(1);
    run@pm.sampling.Sampler(self, getLogFunc, ndim);
end