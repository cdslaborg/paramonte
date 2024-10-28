%>  \brief
%>  Run the ParaDRAM sampler and return nothing.
%>
%>  \details
%>  For example usage, see the documentation of the parent class of this method [pm.sampler.Paradram](@ref Paradram).<br>
%>
%>  \param[inout]   self        :   The input/output parent object of class [pm.sampling.Paradram](@ref Paradram)
%>                                  which is **implicitly** passed to this dynamic method (not by the user).<br>
%>  \param[in]      getLogFunc  :   The input MATLAB function handle or anonymous (lambda) function
%>                                  containing the implementation of the objective function to be sampled.<br>
%>                                  This user-specified function must have the following interface,
%>                                  \code{.m}
%>                                      function logFunc = getLogFunc(state)
%>                                      end
%>                                  \endcode
%>                                  where,
%>                                  <ol>
%>                                      <li>    the input argument ``state`` is a vector of type MATLAB ``double``
%>                                              of size ``ndim`` representing a single point from within the ``ndim``
%>                                              dimensional domain of the mathematical object function to be explored.<br>
%>                                      <li>    the output argument `logFunc` is a scalar of the same type as the
%>                                              input ``state`` containing the natural logarithm of the objective
%>                                              function at the specified input ``state`` within its domain.<br>
%>                                  </ol>
%>  \param[in]      ndim        :   The input scalar positive-valued whole-number representing the number of dimensions
%>                                  of the domain of the user-specified objective function in the input ``getLogFunc()``.
%>
%>  \interface{run}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      sampler.run(getLogFunc, ndim);
%>
%>  \endcode
%>
%>  \final{run}
%>
%>  \author
%>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin%>
function run(self, getLogFunc, ndim)
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
    self.name = string(inputname(1));
    run@pm.sampling.Sampler(self, getLogFunc, ndim);
end