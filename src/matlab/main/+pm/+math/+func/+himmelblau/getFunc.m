%>  \brief
%>  Return the value of the 2-dimensional Himmelblau
%>  function at the specified input value.
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
%>      H(3.0,2.0) = 0.0 ~,<br>
%>      H(-2.805118, 3.131312)=0.0 ~,<br>
%>      H(-3.779310, -3.283186)=0.0 ~,<br>
%>      H(3.584428, -1.848126)=0.0 ~.<br>
%>  \f}
%>
%>  The function is named after David Mautner Himmelblau (1924–2011), who introduced it.
%>  The locations of all the minima can be found analytically.<br>
%>
%>  \param[in]  x   :   The input scalar or array of the same rank and shape as
%>                      other input array-like arguments of type MATLAB ``double``,
%>                      representing the x-component of the state within the domain
%>                      of Himmelblau function at which the function value must be computed.
%>  \param[in]  y   :   The input scalar or array of the same rank and shape as
%>                      other input array-like arguments of type MATLAB ``double``,
%>                      representing the y-component of the state within the domain
%>                      of Himmelblau function at which the function value must be computed.
%>
%>  \return
%>  `func`          :   The output scalar or array of the same rank and shape
%>                      as the input array-like arguments representing the value
%>                      of the Himmelblau function at the specified input ``x`` and ``y``.
%>
%>  \interface{getFunc}
%>  \code{.m}
%>
%>      func = pm.stats.himmelblau.getFunc(x, y);
%>
%>  \endcode
%>
%>  \see
%>  [pm.stats.himmelblau.getFunc](@ref getFunc)<br>
%>  [pm.stats.himmelblau.getLogUDF](@ref getLogUDF)<br>
%>  [pm.sampling.Paradram](@ref Paradram)<br>
%>
%>  \example{getFunc}
%>  \include{lineno} example/stats/himmelblau/getFunc/main.m
%>  \vis{getFunc}
%>  \image html example/stats/himmelblau/getFunc/getFunc.2d.png width=700
%>  \image html example/stats/himmelblau/getFunc/getFunc.3d.png width=700
%>
%>  \final{getFunc}
%>
%>  \author
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin
function func = getFunc(x, y)
    func = (x.^2 + y - 11).^2 + (x + y.^2 - 7).^2;
end
