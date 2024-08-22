%>  \brief
%>  Generate and return the **probabilistically-rounded** to the nearest integer value of the input real number.
%>
%>  \details
%>  Unlike the intrinsic Fortran `nint()` function, the procedures of this generic
%>  interface round the input real number to the nearest integer value **probabilistically**.<br>
%>  The closer the input `val` is to its `nint(val)` result, the more likely it will be rounded to it.<br>
%>  See below for example usage.<br>
%>
%>  \param[in]  val     :   The input scalar or array of the same rank as other array-like arguments,
%>                          containing the number to be rounded probabilistically.<br>
%>
%>  \return
%>  `whole`             :   The output scalar or array of the same rank and shape as the input arguments,
%>                          of default `integer` kind, containing the **probabilistically-rounded** input value to the nearest integer.<br>
%>
%>  \interface{pnint}
%>  \code{.m}
%>
%>      whole = pm.math.pnint(val);
%>
%>  \endcode
%>
%>  \impure
%>
%>  \elemental
%>
%>  \example{pnint}
%>  \include{lineno} example/math/pnint/main.m
%>  \output{pnint}
%>  \include{lineno} example/math/pnint/main.out.m
%>
%>  \final{pnint}
%>
%>  \author
%>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
function whole = pnint(val)
    whole = round(val);
    residual = val - whole; % abs(residual) < .5 always holds.
    urand = unifrnd(0, 1, size(val, 1), size(val, 2));
    mask = 1 - abs(residual) < urand;
    mask1 = mask & residual < 0;
    mask2 = mask & residual >= 0;
    whole(mask1(:)) = whole(mask1(:)) - 1;
    whole(mask2(:)) = whole(mask2(:)) + 1;
end