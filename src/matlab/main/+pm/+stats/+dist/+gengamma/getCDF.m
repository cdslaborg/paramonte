%>  \brief
%>  Return the corresponding Cumulative Distribution Function (CDF)
%>  of the input random value(s) from the Generalized Gamma distribution.<br>
%>
%>  \details
%>  A variable \f$X\f$ is said to be **Generalized Gamma** (GenGamma) distributed if its PDF with the **scale** \f$0 < \sigma < +\infty\f$,
%>  **shape** \f$0 < \omega < +\infty\f$, and **shape** \f$0 < \kappa < +\infty\f$ parameters is described by the following equation,
%>
%>  \f{equation}{
%>      \large
%>      \pi(x | \kappa, \omega, \sigma) =
%>      \frac{1}{\sigma \omega \Gamma(\kappa)} ~
%>      \bigg( \frac{x}{\sigma} \bigg)^{\frac{\kappa}{\omega} - 1} \exp\Bigg( -\bigg(\frac{x}{\sigma}\bigg)^{\frac{1}{\omega}} \Bigg) ~,~ 0 < x < \infty
%>  \f}
%>
%>  where \f$\eta = \frac{1}{\sigma \omega \Gamma(\kappa)}\f$ is the **normalization factor** of the PDF.<br>
%>  When \f$\sigma = 1\f$, the GenGamma PDF simplifies to the form,
%>
%>  \f{equation}{
%>      \large
%>      \pi(x) =
%>      \frac{1}{\omega \Gamma(\kappa)} ~
%>      x^{\frac{\kappa}{\omega} - 1} \exp\Bigg( -x^{\frac{1}{\omega}} \Bigg) ~,~ 0 < x < \infty
%>  \f}
%>
%>  If \f$(\sigma, \omega) = (1, 1)\f$, the GenGamma PDF further simplifies to the form,
%>
%>  \f{equation}{
%>      \large
%>      \pi(x) =
%>      \frac{1}{\Gamma(\kappa)} ~
%>      x^{\kappa - 1} \exp(-x) ~,~ 0 < x < \infty
%>  \f}
%>
%>  Setting the shape parameter to \f$\kappa = 1\f$ further simplifies the PDF to the Exponential distribution PDF with the scale parameter \f$\sigma = 1\f$,
%>
%>  \f{equation}{
%>      \large
%>      \pi(x) = \exp(x) ~,~ 0 < x < \infty
%>  \f}
%>
%>  <ol>
%>      <li>    The parameter \f$\sigma\f$ determines the scale of the GenGamma PDF.<br>
%>      <li>    When \f$\omega = 1\f$, the GenGamma PDF reduces to the PDF of the [Gamma distribution](@ref pm_distGamma).<br>
%>      <li>    When \f$\kappa = 1, \omega = 1\f$, the GenGamma PDF reduces to the PDF of the [Exponential distribution](@ref pm_distExp).<br>
%>  </ol>
%>
%>  The **CDF** of the Generalized Gamma distribution over a strictly-positive support \f$x \in (0, +\infty)\f$ with the three (shape, shape, scale)
%>  parameters \f$(\kappa > 0, \omega > 0, \sigma > 0)\f$ is defined by the **regularized** [Lower Incomplete Gamma function](@ref pm_mathGamma) as,
%>  \f{eqnarray}{
%>      \large
%>      \mathrm{CDF}(x | \kappa, \omega, \sigma)
%>      & = & P\bigg(\kappa, \big(\frac{x}{\sigma}\big)^{\frac{1}{\omega}} \bigg) \\
%>      & = & \frac{1}{\Gamma(\kappa)} \int_0^{\big(\frac{x}{\sigma}\big)^{\frac{1}{\omega}}} ~ t^{\kappa - 1}{\mathrm e}^{-t} ~ dt ~,
%>  \f}
%>
%>  where \f$\Gamma(\kappa)\f$ represents the Gamma function.<br>
%>
%>  \param[in]  x           :   The input scalar or array of the same shape as other array-valued arguments,
%>                              containing the values at which the CDF must be computed.<br>
%>  \param[in]  kappa       :   The input scalar or array of the same shape as other array-valued arguments,
%>                              containing the shape parameter (\f$\kappa\f$) of the distribution.<br>
%>  \param[in]  invOmega    :   The input scalar or array of the same shape as other array-valued arguments,
%>                              containing the inverse of the second shape parameter (\f$\omega\f$) of the distribution.<br>
%>  \param[in]  invSigma    :   The input scalar or array of the same shape as other array-valued arguments,
%>                              containing the inverse of the scale parameter (\f$\sigma\f$) of the distribution.<br>
%>
%>  \return
%>  `cdf`                :   The output scalar or array of the same shape as other array-valued input arguments
%>                              of the same type and kind as `x` containing the PDF of the specified GenGamma distribution.
%>
%>  \interface{getCDF}
%>  \code{.m}
%>
%>      cdf = pm.stats.dist.gengamma.getCDF(x, kappa, invOmega, invSigma)
%>
%>  \endcode
%>
%>  \example{getCDF}
%>  \include{lineno} example/stats/dist/gengamma/getCDF/main.m
%>  \output{getCDF}
%>  \include{lineno} example/stats/dist/gengamma/getCDF/main.out.m
%>  \vis{getCDF}
%>  \image html example/stats/dist/gengamma/getCDF/gengamma.getCDF.line.png width=700
%>
%>  \final{getCDF}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, July 19, 2024, 1:07 AM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function cdf = getCDF(x, kappa, invOmega, invSigma)
    assert(all(0 <= x), 'The condition `all(0 < x)` must hold.');
    assert(all(0 < kappa), 'The condition `all(0 < kappa)` must hold.');
    assert(all(0 < invOmega), 'The condition `all(0 < invOmega)` must hold.');
    assert(all(0 < invSigma), 'The condition `all(0 < invSigma)` must hold.');
    cdf = gammainc((x .* invSigma).^invOmega, kappa);
end