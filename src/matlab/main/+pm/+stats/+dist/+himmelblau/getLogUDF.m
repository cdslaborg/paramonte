%>  \brief
%>  Return the natural logarithm of the Unnormalized Density Function (UDF)
%>  of the inverse of the 2-dimensional modified Himmelblau function.
%>
%>  \details
%>  Himmelblau's function is a multi-modal function, used to test
%>  the performance of optimization algorithms. The function is defined by:
%>  \f{equation}
%>      H(x, y) = (x^{2} + y - 11)^{2} + (x + y^{2} - 7)^{2} ~.
%>  \f}
%>
%>  It has one local maximum at \f$x = -0.270845\f$ and \f$y = -0.923039\f$ where \f$H(x, y) = 181.617\f$,
%>  and four identical local minima:
%>  \f{equation}
%>      H(3.0,2.0) = 0.0 ~,
%>      H(-2.805118, 3.131312)=0.0 ~,
%>      H(-3.779310, -3.283186)=0.0 ~,
%>      H(3.584428, -1.848126)=0.0 ~.
%>  \f}
%>
%>  The function is named after David Mautner Himmelblau (1924–2011), who introduced it.<br>
%>  The locations of all the minima can be found analytically.<br>
%>
%>  This MATLAB function returns a modification of the Himmelblau function as a density
%>  function suitable for testing sampling algorithms (or stochastic maximizers):
%>  \f{equation}
%>      f(x, y, \epsilon) = \frac{1}{H(x, y) + \epsilon} ~.
%>  \f}
%>
%>  \param[in]  x       :   The input scalar or array of the same rank and shape as
%>                          other input array-like arguments of type MATLAB ``double``,
%>                          representing the x-component of the state within the domain
%>                          of Himmelblau density at which the density value must be computed.<br>
%>  \param[in]  y       :   The input scalar or array of the same rank and shape as
%>                          other input array-like arguments of type MATLAB ``double``,
%>                          representing the y-component of the state within the domain
%>                          of Himmelblau density at which the density value must be computed.<br>
%>
%>  \param[in]  epsilon :   The input positive scalar or array of the same rank and shape as
%>                          other array-like input arguments of type MATLAB ``double``,
%>                          representing the value to be added to the inverse
%>                          of the Himmelblau function.<br>
%>                          Increasingly smaller values of ``epsilon`` will yield pointier densities.<br>
%>                          Increasingly larger values of ``epsilon`` will yield flatter densities.<br>
%>                          (**optional**, default = ``1``)
%>
%>  \return
%>  `logUDF`            :   The output Unnormalized Density Function (UDF) of the
%>                          Himmelblau density at the specified input ``state``.
%>
%>  \interface{getLogUDF}
%>  \code{.m}
%>
%>      logUDF = pm.stats.himmelblau.getLogUDF(state);
%>      logUDF = pm.stats.himmelblau.getLogUDF(state, epsilon);
%>
%>  \endcode
%>
%>  \see
%>  [pm.stats.himmelblau.getFunc](@ref getFunc)<br>
%>  [pm.stats.himmelblau.getLogUDF](@ref getLogUDF)<br>
%>  [pm.sampling.Paradram](@ref Paradram)<br>
%>
%>  \example{getLogUDF}
%>  \include{lineno} example/stats/himmelblau/getLogUDF/main.m
%>  \vis{getLogUDF}
%>  \image html example/stats/himmelblau/getLogUDF/getLogUDF.2d.png width=700
%>  \image html example/stats/himmelblau/getLogUDF/getLogUDF.3d.png width=700
%>
%>  \final{getLogUDF}
%>
%>  \author
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin
function logFunc = getLogUDF(x, y, epsilon)
    if  nargin < 3
        epsilon = 1;
    end
    logFunc = -log(pm.math.func.himmelblau.getFunc(x, y) + epsilon);
    %if false
    %    nsim = 100000;
    %    for i = 1 : nsim
    %        logFunc = logFunc - log((state(1)^2 + state(2) - 11)^2 + (state(1) + state(2)^2 - 7)^2 + 0.1);
    %    end
    %    logFunc = logFunc / (nsim + 1);
    %end
end