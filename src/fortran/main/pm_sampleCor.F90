!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains classes and procedures for computing properties related to the correlation matrices of random samples.
!>
!>  \details
!>
!>  Correlation matrix
!>  ==================
!>
!>  The correlation matrix of \f$N\f$ random variables \f$X_{1}, \ldots, X_{N}\f$
!>  is the \f$N\times N\f$ matrix \f$\rho\f$ whose \f$(i, j)\f$ entry is,
!>  \f{equation}{
!>      \rho_{ij} := \up{COR}(X_{i}, X_{j}) = \frac{\up{COV}(X_{i}, X_{j})}{\sigma_{X_{i}} \sigma_{X_{j}}}, \quad {\text{if}} ~ \sigma_{X_{i}} \sigma_{X_{j}} > 0 ~.
!>  \f}
!>
!>  Thus the diagonal entries are all identically one.<br>
!>  If the measures of correlation used are product-moment coefficients, the correlation matrix is the same as the covariance matrix of the standardized random variables
!>  \f$X_{i} / \sigma(X_{i})\f$ for \f$i = 1, \dots, N\f$.<br>
!>  This applies both to the matrix of population correlations (in which case \f$\sigma\f$ is the population standard deviation),
!>  and to the matrix of sample correlations (in which case \f$\sigma\f$ denotes the sample standard deviation).<br>
!>  Consequently, each is necessarily a positive-semidefinite matrix.<br>
!>  Moreover, the correlation matrix is strictly positive definite if no variable
!>  can have all its values exactly generated as a linear function of the values of the others.<br>
!>  The correlation matrix is symmetric because the correlation between \f$X_{i}\f$ and \f$X_{j}\f$ is the same as the correlation between \f$X_{j}\f$ and \f$X_{i}\f$.<br>
!>
!>  \note
!>  The best way to compute the correlation matrix of a sample is to first compute the covariance matrix of the sample via the
!>  the relevant procedures in [pm_sampleCov](@ref pm_sampleCov) and then call the relevant procedures of this module to
!>  convert the covariance matrix to the corresponding correlation matrix.<br>
!>  The elements of the correlation matrix can be computed via the following equation,
!>  \f{equation}{
!>      \rho_{ij} = \frac{\Sigma_{ij}} { \sqrt{\Sigma_{ii}} \times \sqrt{\Sigma_{jj}} } ~,
!>  \f}
!>  where \f$\rho\f$ represents the correlation matrix, \f$\Sigma\f$ represents the covariance matrix, and \f$\sigma\f$ represents the standard deviation.<br>
!>
!>  \note
!>  The vector of standard deviations can be readily extracted from the covariance matrix via [getMatCopy](@ref pm_matrixCopy::getMatCopy).<br>
!>
!>  Sample Cordance
!>  ===============
!>
!>  Let \f$\{(X_{1}, Y_{1}), \ldots, (X_{N}, Y_{N})\}\f$ be a set of \f$N\f$ observations of the joint random variables \f$X\f$ and \f$Y\f$.<br>
!>  A concordant pair is a pair of observations, each on two variables, \f$(X_i, Y_i)\f$ and \f$(X_j, Y2_j)\f$, having the property that
!>  \f{equation}{
!>      \up{sgn}(X_j - X_i) ~=~ \up{sgn}(Y_j - Y_i) ~,
!>  \f}
!>  where \f$\up{sgn}\f$ is the [signum function](https://en.wikipedia.org/wiki/Sign_function) defined as:<br>
!>  \f{equation}{
!>      \up{sgn}(x) =
!>      \begin{cases}
!>          -1, & x < 0 ~, \\
!>           0, & x = 0 ~, \\
!>           1, & x > 0 ~,
!>      \end{cases}
!>  \f}
!>  that is, in a concordant pair, both elements of one pair are either greater than, equal to, or less than the corresponding elements of the other pair.<br>
!>  In contrast, a **discordant pair** is a pair of two-variable observations such that,
!>  \f{equation}{
!>      \up{sgn} (X_j - X_i) ~=~ -\up{sgn}(Y_j - Y_i) ~,
!>  \f}
!>  that is, if one pair contains a higher value of \f$X\f$ then the **other** pair contains a higher value of \f$Y\f$.<br>
!>
!>  Sample concordance is relevant to computing the [Kendall correlation coefficient](@ref pm_sampleCor) or in hypothesis testing.<br>
!>  However, in many situations it is also important to distinguish **tied pairs** from concordant pairs.<br>
!>  Therefore, more precisely, any pair of observations \f$(X_{i}, Y_{i})\f$ and \f$(X_{j}, Y_{j})\f$ are considered,
!>  <ol>
!>      <li>    **X-tied** if \f$\up{sgn} (X_j - X_i) = 0\f$.<br>
!>      <li>    **Y-tied** if \f$\up{sgn} (Y_j - Y_i) = 0\f$.<br>
!>      <li>    **concordant** if \f$\ms{sgn}(X_{i} - X_{j}) = \ms{sgn}(Y_{i} - Y_{j})\f$
!>              where \f$\ms{sgn}(\cdots)\f$ is the [sign function](https://en.wikipedia.org/wiki/Sign_function).<br>
!>      <li>    **discordant** if \f$\ms{sgn}(X_{i} - X_{j}) \neq \ms{sgn}(Y_{i} - Y_{j})\f$.<br>
!>  </ol>
!>
!>  The generic interface [setCordance](@ref pm_sampleCor::setCordance) of this module returns 
!>  the **sample cordance tuple/vector** comprised of the number of x-ties, y-ties, concordant pairs, and discordant pairs.<br>
!>
!>  Sample Cordance Algorithm
!>  -------------------------
!>
!>  The naive method of computing sample concordance quickly becomes expensive for large \f$N\f$,
!>  because there are \f${N \choose 2} = \frac{N (N - 1)}{2}\f$ (the binomial coefficient) number of ways to choose two items from \f$N\f$ items.<br>
!>  This makes the naive algorithm of complexity \f$\mathcal{O}(N^{2})\f$.<br>
!>
!>  Correlation coefficient
!>  =======================
!>
!>  A correlation coefficient is a numerical measure of some type of correlation, meaning a statistical relationship between two variables.<br>
!>  The variables may be two columns of a given data set of observations, often called a sample, or two components of a multivariate random variable with a known distribution.<br>
!>
!>  This module contains generic algorithms for computing the following popular sample correlation coefficients.<br>
!>
!>  Pearson correlation coefficient
!>  -------------------------------
!>
!>  The Pearson product-moment correlation coefficient, also known as \f$r\f$, \f$R\f$, or the Pearson \f$r\f$,
!>  is a measure of the strength and direction of the linear relationship between two variables that is defined as
!>  the covariance of the variables divided by the product of their standard deviations.<br>
!>  This is the best-known and most commonly used type of **correlation coefficient**.<br>
!>  When the term *correlation coefficient* is used without further qualification,
!>  it usually refers to the **Pearson product-moment correlation coefficient**.<br>
!>
!>  **Definition**<br>
!>
!>  The Pearson correlation coefficient, when applied to a **population**, is commonly represented by the Greek letter \f$\rho\f$ (rho) and
!>  may be referred to as the population correlation coefficient or the population Pearson correlation coefficient.<br>
!>  Given a pair of random variables \f$(X,Y)\f$, the formula for \f$\rho\f$ is,<br>
!>  \f{equation}{
!>      \rho_{X,Y}={\frac {\up{cov} (X,Y)}{\sigma_{X}\sigma_{Y}}} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$\up{cov}\f$ is the covariance.
!>      <li>    \f$\sigma_{X}\f$ is the standard deviation of \f$X\f$.
!>      <li>    \f$\sigma_Y\f$ is the standard deviation of \f$Y\f$.
!>  </ol>
!>
!>  The formula for \f$\rho\f$ can be also expressed in terms of mean and expectation. Since,
!>  \f{equation}{
!>      \up{cov}(X,Y) = \up{\mathbb {E}} [(X-\mu _{X})(Y-\mu _{Y})] ~,
!>  \f}
!>  the formula for \f$\rho\f$ can also be written as,
!>  \f{equation}{
!>      \rho_{X,Y} = {\frac{\up{\mathbb{E}} [(X - \mu_{X})(Y - \mu_{Y})]}{\sigma _{X}\sigma _{Y}}} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$\sigma_Y\f$ and \f$\sigma_{X}\f$ are defined as above.
!>      <li>    \f$\mu_{X}\f$ is the mean of \f$X\f$.
!>      <li>    \f$\mu_{Y}\f$ is the mean of \f$Y\f$.
!>      <li>    \f$\up{\mathbb{E}}\f$ is the expectation.
!>  </ol>
!>
!>  The Pearson correlation coefficient, when applied to a **sample**, is commonly represented by \f$r_{xy}\f$ and
!>  may be referred to as the **sample correlation coefficient** or the **sample Pearson correlation coefficient**.<br>
!>  We can obtain a formula for \f$r_{xy}\f$ by substituting estimates of the covariances and variances based on a sample into the formula above.<br>
!>  Given paired data \f$\left\{(x_{1}, y_{1}), \ldots,(x_{n}, y_{n})\right\}\f$ consisting of \f$n\f$ pairs, \f$r_{xy}\f$ is defined as,<br>
!>  \f{equation}{
!>      r_{xy} = {\frac {\sum_{i=1}^{n}(x_{i}-{\bar {x}})(y_{i}-{\bar {y}})}{{\sqrt {\sum _{i=1}^{n}(x_{i}-{\bar {x}})^{2}}}{\sqrt {\sum _{i=1}^{n}(y_{i}-{\bar {y}})^{2}}}}} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$n\f$ is sample size.
!>      <li>    \f$x_{i},y_{i}\f$ are the individual sample points indexed with \f$i\f$.
!>      <li>    \f$\textstyle {\bar{x}} = {\frac{1}{n}} \sum_{i=1}^{n}x_{i}\f$ (the sample mean); and analogously for \f$\bar{y}\f$.
!>  </ol>
!>
!>  Spearman rank correlation coefficient
!>  -------------------------------------
!>
!>  The Spearman correlation coefficient is defined as the Pearson correlation coefficient between the **rank variables**.<br>
!>  For a sample of size \f$n\f$, the \f$n\f$ raw scores \f$X_{i}, Y_{i}\f$ are converted to ranks \f$\up{R}({X_{i}}), \up{R}({Y_{i}})\f$, and \f$r_{s}\f$ is computed as
!>  \f{equation}{
!>      r_{s} = \rho_{\up {R} (X),\up {R} (Y)} = {\frac {\up {cov} (\up {R} (X),\up {R} (Y))}{\sigma _{\up {R} (X)}\sigma _{\up {R} (Y)}}} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$\rho\f$ denotes the usual Pearson correlation coefficient, but applied to the rank variables,
!>      <li>    \f$\up{cov} (\up {R} (X),\up {R} (Y))\f$ is the covariance of the rank variables,
!>      <li>    \f$\sigma_{\up {R} (X)}\f$ and \f$\sigma_{\up {R} (Y)}\f$ are the standard deviations of the rank variables.
!>  </ol>
!>
!>  When all \f$n\f$ ranks are **distinct integers**, it can be computed using the formula,<br>
!>  \f{equation}{
!>      r_{s} = 1 - {\frac{6\sum d_{i}^{2}}{n(n^{2}-1)}} ~,
!>  \f}
!>  where
!>  <ol>
!>      <li>    \f$d_{i}=\up {R} (X_{i})-\up{R} (Y_{i})\f$ is the difference between the two ranks of each observation,
!>      <li>    \f$n\f$ is the number of observations.
!>  </ol>
!>
!>  Kendall rank correlation coefficient
!>  ------------------------------------
!>
!>  The Kendall rank correlation coefficient, commonly referred to as the Kendall \f$\tau\f$ coefficient,
!>  is a statistic used to measure the ordinal association between two measured quantities.<br>
!>  It is a measure of **rank correlation**: the similarity of the orderings of the data when ranked by each of the quantities.<br>
!>  It is named after Maurice Kendall, who developed it in 1938, though Gustav Fechner had proposed a similar measure in the context of time series in 1897.<br>
!>
!>  **Definition**<br>
!>
!>  Let \f$(x_{1},y_{1}), \ldots,(x_{n},y_{n})\f$ be a set of observations of the joint random variables \f$X\f$ and \f$Y\f$,
!>  such that all the values \f$x_{i}\f$ and \f$y_{i}\f$ are unique.<br>
!>  Any pair of observations \f$(x_{i},y_{i})\f$ and \f$(x_{j},y_{j})\f$, where \f$i < j\f$, are said to be **concordant** if the sort order of \f$(x_{i}, x_{j})\f$ and \f$(y_{i},y_{j})\f$ agrees,<br>
!>  that is, if either both \f$x_{i} > x_{j}\f$ and \f$y_{i} > y_{j}\f$ holds or both \f$x_{i} < x_{j}\f$ and \f$y_{i} < y_{j}\f$; otherwise they are said to be **discordant**.<br>
!>  The Kendall \f$\tau\f$ coefficient is defined as:<br>
!>  \f{equation}{
!>      \tau = {\frac {({\text{number of concordant pairs}})-({\text{number of discordant pairs}})}{({\text{number of pairs}})}} = 1 - {\frac {2({\text{number of discordant pairs}})}{n \choose 2}} ~.
!>  \f}
!>  where \f${n \choose 2} = \frac{n(n-1)}{2}\f$ is the binomial coefficient for the number of ways to choose two items from \f$n\f$ items.<br>
!>
!>  **Accounting for ties**<br>
!>
!>  A pair \f$\{(x_{i}, y_{i}), (x_{j}, y_{j})\}\f$ is said to be tied if and only if \f$x_{i} = x_{j}\f$ or \f$y_{i} = y_{j}\f$.<br>
!>  A tied pair is neither concordant nor discordant.<br>
!>  When tied pairs arise in the data, the coefficient may be modified in a number of ways to keep it in the range \f$[−1, 1]\f$:<br>
!>
!>  **Tau-a**<br>
!>
!>  The Tau-a statistic tests the strength of association of the cross tabulations.<br>
!>  Both variables have to be ordinal.<br>
!>  The Tau-a coefficient does **not** make any adjustment for ties.<br>
!>  It is defined as:<br>
!>  \f{equation}{
!>      \tau_{A} = {\frac {n_{c}-n_{d}}{n_{0}}} ~,
!>  \f}
!>  where \f$n_c\f$, \f$n_d\f$ and \f$n_0\f$ are defined as in below for Tau-b coefficient.<br>
!>
!>  **Tau-b**<br>
!>
!>  The Tau-b statistic, unlike Tau-a, makes adjustments for ties.<br>
!>  Values of Tau-b range from \f$−1\f$ (\f$100\%\f$ negative association, or perfect inversion) to \f$+1\f$ (\f$100\%\f$ positive association, or perfect agreement).<br>
!>  A value of zero indicates the absence of association.<br>
!>  The Kendall Tau-b coefficient is defined as:<br>
!>  \f{equation}{
!>      \tau_{B} = \frac{n_{c}-n_{d}}{\sqrt {(n_{0}-n_{1})(n_{0}-n_{2})}} ~,
!>  \f}
!>  where
!>  \f{aligned}{
!>      n_{0} &= n(n-1)/2 \\
!>      n_{1} &= \sum _{i}t_{i}(t_{i}-1)/2 \\
!>      n_{2} &= \sum _{j}u_{j}(u_{j}-1)/2 \\
!>      n_{c} &= \text{Number of concordant pairs} \\
!>      n_{d} &= \text{Number of discordant pairs} \\
!>      t_{i} &= \text{Number of tied values in the } i^{\text{th}} \text{ group of ties for the first quantity} \\
!>      u_{j} &= \text{Number of tied values in the } j^{\text{th}} \text{ group of ties for the second quantity}
!>  \f}
!>
!>  There are also other definitions of Kendall \f$\tau\f$ which are not considered in the current version of the ParaMonte library.<br>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>  [pm_sampleACT](@ref pm_sampleACT)<br>
!>  [pm_sampleCCF](@ref pm_sampleCCF)<br>
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>  [pm_sampleCov](@ref pm_sampleCov)<br>
!>  [pm_sampleConv](@ref pm_sampleConv)<br>
!>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
!>  [pm_sampleMean](@ref pm_sampleMean)<br>
!>  [pm_sampleNorm](@ref pm_sampleNorm)<br>
!>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
!>  [pm_sampleScale](@ref pm_sampleScale)<br>
!>  [pm_sampleShift](@ref pm_sampleShift)<br>
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_sampleAffinity](@ref pm_sampleAffinity)<br>
!>  [pm_sampleVar](@ref pm_sampleVar)<br>
!>  [Correlation coefficient](https://en.wikipedia.org/wiki/Correlation_coefficient)<br>
!>
!>  \test
!>  [test_pm_sampleCor](@ref test_pm_sampleCor)
!>
!>  \todo
!>  \phigh
!>  The procedures of this module should be extended to support samples of type `complex` of arbitrary kind type parameter, similar to procedures of [pm_sampleVar](@ref pm_sampleVar).<br>
!>
!>  \final{pm_sampleCor}
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 01:45 AM, August 21, 2018, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleCor

    use pm_kind, only: SK, IK, LK, RK
    use pm_matrixSubset, only: uppDia_type, uppDia
    use pm_matrixSubset, only: lowDia_type, lowDia
    use pm_matrixSubset, only: upp_type, upp
    use pm_matrixSubset, only: low_type, low
    use pm_container, only: css_type
#if PDT_ENABLED
    use pm_container, only: css_pdt
#endif

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleCor"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different correlation coefficients (e.g., pearson, spearman, kendall, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{corcoef_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract :: corcoef_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the **kendall** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [kendall](@ref pm_sampleCor::kendall) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{kendall_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(corcoef_type) :: kendall_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [kendall_type](@ref pm_sampleCor::kendall_type) that is exclusively used
    !>  to signify the **kendall** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{kendall}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(kendall_type), parameter :: kendall = kendall_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: kendall
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the **kendallA** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [kendallA](@ref pm_sampleCor::kendallA) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [kendallA](@ref pm_sampleCor::kendallA)<br>
    !>  [kendallB](@ref pm_sampleCor::kendallB)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [kendallA_type](@ref pm_sampleCor::kendallA_type)<br>
    !>  [kendallB_type](@ref pm_sampleCor::kendallB_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{kendallA_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(kendall_type) :: kendallA_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [kendallA_type](@ref pm_sampleCor::kendallA_type) that is exclusively used
    !>  to signify the **kendallA** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [kendallA](@ref pm_sampleCor::kendallA)<br>
    !>  [kendallB](@ref pm_sampleCor::kendallB)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [kendallA_type](@ref pm_sampleCor::kendallA_type)<br>
    !>  [kendallB_type](@ref pm_sampleCor::kendallB_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{kendallA}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(kendallA_type), parameter :: kendallA = kendallA_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: kendallA
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the **kendallB** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [kendallB](@ref pm_sampleCor::kendallB) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [kendallA](@ref pm_sampleCor::kendallA)<br>
    !>  [kendallB](@ref pm_sampleCor::kendallB)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [kendallA_type](@ref pm_sampleCor::kendallA_type)<br>
    !>  [kendallB_type](@ref pm_sampleCor::kendallB_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{kendallB_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(kendall_type) :: kendallB_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [kendallB_type](@ref pm_sampleCor::kendallB_type) that is exclusively used
    !>  to signify the **kendallB** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [kendallA](@ref pm_sampleCor::kendallA)<br>
    !>  [kendallB](@ref pm_sampleCor::kendallB)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [kendallA_type](@ref pm_sampleCor::kendallA_type)<br>
    !>  [kendallB_type](@ref pm_sampleCor::kendallB_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{kendallB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(kendallB_type), parameter :: kendallB = kendallB_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: kendallB
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the **pearson** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  Frequency weights, represent the number of duplications of each observation in the sample.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [pearson](@ref pm_sampleCor::pearson) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{pearson_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(corcoef_type) :: pearson_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [pearson_type](@ref pm_sampleCor::pearson_type) that is exclusively used
    !>  to signify the **pearson** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{pearson}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(pearson_type), parameter :: pearson = pearson_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: pearson
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the **spearman** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [spearman](@ref pm_sampleCor::spearman) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{spearman_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(corcoef_type) :: spearman_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [spearman_type](@ref pm_sampleCor::spearman_type) that is exclusively used
    !>  to signify the **spearman** type of correlation coefficients.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [kendall](@ref pm_sampleCor::kendall)<br>
    !>  [pearson](@ref pm_sampleCor::pearson)<br>
    !>  [spearman](@ref pm_sampleCor::spearman)<br>
    !>  [corcoef_type](@ref pm_sampleCor::corcoef_type)<br>
    !>  [kendall_type](@ref pm_sampleCor::kendall_type)<br>
    !>  [pearson_type](@ref pm_sampleCor::pearson_type)<br>
    !>  [spearman_type](@ref pm_sampleCor::spearman_type)<br>
    !>
    !>  \final{spearman}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(spearman_type), parameter :: spearman = spearman_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: spearman
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (Pearson) correlation coefficient or matrix of a pair of (weighted) time series `x(1:nsam)` and `y(1:nsam)`
    !>  or of an input (weighted) array of shape `(ndim, nsam)` or `(nsam, ndim)` where `ndim` is the number of data dimensions
    !>  (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  This generic interface performs one of the following computational tasks:<br>
    !>  <ol>
    !>      <li>    compute the correlation matrix corresponding to an input covariance matrix and vector of standard deviations in arbitrary `ndim` dimensions.<br>
    !>      <li>    compute the sample correlation matrix of a random sample of `nsam` observations in arbitrary `ndim` dimensional space.<br>
    !>      <li>    compute the sample correlation matrix of a pair of `nsam` observations in two separated data vectors `x` and `y`.<br>
    !>  </ol>
    !>  See the documentation of the parent module [pm_sampleCor](@ref pm_sampleCor) for algorithmic details and sample correlation matrix definition.<br>
    !>
    !>  \note
    !>  The **sample correlation matrix** can be also readily computed from the [sample covariance matrix](@ref pm_sampleCov) as in the following example,
    !>  \code{.F90}
    !>      cor(:,:) = getCor(getCov(sample, dim))
    !>      cor(:,:) = getCor(getCov(sample, dim, weight)) ! weighted sample.
    !>      !
    !>  \endcode
    !>
    !>  \param[in]  cov         :   The input positive semi-definite square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the output `cor`,
    !>                              representing the correlation matrix based upon which the output correlation matrix `cor` is to be computed.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input arguments `sample` and `x` and `y` are missing.)
    !>  \param[in]  subsetv     :   The input scalar constant argument that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input covariance matrix must be used.<br>
    !>                                  <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input covariance matrix must be used.<br>
    !>                                  <li>    The constant [low](@ref pm_matrixSubset::low), implying that only the lower subset of the input covariance matrix must be used.<br>
    !>                                          This option is available **if and only if** the input argument `stdinv` is present (otherwise, how do we know the variances?).<br>
    !>                                  <li>    The constant [upp](@ref pm_matrixSubset::upp), implying that only the upper subset of the input covariance matrix must be used.<br>
    !>                                          This option is available **if and only if** the input argument `stdinv` is present (otherwise, how do we know the variances?).<br>
    !>                              </ol>
    !>                              This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `cov` is present.)
    !>  \param[in]  stdinv      :   The input positive vector of shape `(1 : ndim)` of type `real` of the same kind as the input `cor`,
    !>                              containing the inverse of the standard deviations of the `ndim` data attributes based upon which the output correlation matrix `cor` is to be computed.<br>
    !>                              \f{equation}{
    !>                                  \ms{stdinv}(i) = \frac{1} {\sqrt{\ms{cov}(i,i)}} ~,
    !>                              \f}
    !>                              (**optional**. It must be present **if and only if** the input argument `subsetv` is set to either [low](@ref pm_matrixSubset::low) or [upp](@ref pm_matrixSubset::upp).)
    !>  \param[in]  x           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cor`,
    !>                              containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `y` is present and `sample` is missing.)
    !>  \param[in]  y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cor`,
    !>                              containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `x` is present and `sample` is missing.)
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `cor`,
    !>                              containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                              If `sample` is a matrix, then the direction along which the correlation matrix is computed is dictated by the input argument `dim`.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `dim` is present and `x` and `y` are missing.)
    !>  \param[in]  dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the correlation matrix must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input argument `sample` is present.)
    !>  \param[in]  weight      :   The `contiguous` vector of length `nsam` of,
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK, or
    !>                                  <li>    type `real` of the same kind as the kind of the output `cor`,
    !>                              </ol>
    !>                              containing the corresponding weights of individual `nsam` observations in `sample`.<br>
    !>                              (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled). It can be present **only if** the input arguments `cov`, `subsetv`, and `stdinv` are missing.)
    !>
    !>  \return
    !>  `cor`                   :   The output **positive semi-definite** square matrix of shape `(1:ndim, 1:ndim)` of,
    !>                              <ol>
    !>                                  <li>     type `complex` of kind \CKALL,
    !>                                  <li>     type `real` of kind \RKALL,
    !>                              </ol>
    !>                              representing its Pearson correlation matrix.<br>
    !>                              When the input arguments `x` and `y` are present, then `ndim == 2` holds by definition and the output `cor` is a scalar of value \f$r_{xy}\f$,
    !>                              \f{equation}{
    !>                                  \ms{cor} =
    !>                                  \begin{bmatrix}
    !>                                      1 && r_{xy} \\
    !>                                      r_{yx} && 1 ~.
    !>                                  \end{bmatrix}
    !>                              \f}
    !>
    !>  \interface{getCor}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCor, only: getCor
    !>
    !>      ! correlation matrix to correlation matrix.
    !>
    !>      cor(1:ndim, 1:ndim) = getCor(cov(1:ndim, 1:ndim), subsetv) ! subsetv = uppDia, lowDia
    !>      cor(1:ndim, 1:ndim) = getCor(cov(1:ndim, 1:ndim), subsetv, stdinv(1:ndim)) ! subsetv = upp, low
    !>
    !>      ! XY time series correlation matrix.
    !>
    !>      cor = getCor(x(1:nsam), y(1:nsam)) ! correlation coefficient.
    !>      cor = getCor(x(1:nsam), y(1:nsam), weight(1:nsam)) ! correlation coefficient.
    !>
    !>      ! sample correlation matrix.
    !>
    !>      cor(:,:) = getCor(sample(:,:), dim) ! full correlation matrix.
    !>      cor(:,:) = getCor(sample(:,:), dim, weight(1:nsam)) ! upper/lower correlation matrix.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  This generic interface is merely a flexible wrapper around the generic `subroutine` interface [setCor](@ref pm_sampleCor::setCor).<br>
    !>  As such, all conditions and warnings associated with [setCor](@ref pm_sampleCor::setCor) equally hold for this generic interface.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{getCor}
    !>  \include{lineno} example/pm_sampleCor/getCor/main.F90
    !>  \compilef{getCor}
    !>  \output{getCor}
    !>  \include{lineno} example/pm_sampleCor/getCor/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
    !>
    !>  \final{getCor}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! cov2cor.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCFC_RULD_VUXD_CK5(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCFC_RULD_VUXD_CK4(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCFC_RULD_VUXD_CK3(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCFC_RULD_VUXD_CK2(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCFC_RULD_VUXD_CK1(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCFC_RULD_VUXD_RK5(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCFC_RULD_VUXD_RK4(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCFC_RULD_VUXD_RK3(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCFC_RULD_VUXD_RK2(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCFC_RULD_VUXD_RK1(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCFC_RULD_VXLD_CK5(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCFC_RULD_VXLD_CK4(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCFC_RULD_VXLD_CK3(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCFC_RULD_VXLD_CK2(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCFC_RULD_VXLD_CK1(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCFC_RULD_VXLD_RK5(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCFC_RULD_VXLD_RK4(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCFC_RULD_VXLD_RK3(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCFC_RULD_VXLD_RK2(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCFC_RULD_VXLD_RK1(cov, subsetv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCFC_RULD_VUXX_CK5(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCFC_RULD_VUXX_CK4(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCFC_RULD_VUXX_CK3(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCFC_RULD_VUXX_CK2(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCFC_RULD_VUXX_CK1(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCFC_RULD_VUXX_RK5(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCFC_RULD_VUXX_RK4(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCFC_RULD_VUXX_RK3(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCFC_RULD_VUXX_RK2(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCFC_RULD_VUXX_RK1(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VUXX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCFC_RULD_VXLX_CK5(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCFC_RULD_VXLX_CK4(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCFC_RULD_VXLX_CK3(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCFC_RULD_VXLX_CK2(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCFC_RULD_VXLX_CK1(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)                                                :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCFC_RULD_VXLX_RK5(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCFC_RULD_VXLX_RK4(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCFC_RULD_VXLX_RK3(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCFC_RULD_VXLX_RK2(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCFC_RULD_VXLX_RK1(cov, subsetv, stdinv) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCFC_RULD_VXLX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)                                                   :: cor(size(cov, 1, IK), size(cov, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! no weight: XY.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPrsWNO_XY_CK5(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK4_ENABLED
    PURE module function getPrsWNO_XY_CK4(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK3_ENABLED
    PURE module function getPrsWNO_XY_CK3(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK2_ENABLED
    PURE module function getPrsWNO_XY_CK2(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK1_ENABLED
    PURE module function getPrsWNO_XY_CK1(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPrsWNO_XY_RK5(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK4_ENABLED
    PURE module function getPrsWNO_XY_RK4(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK3_ENABLED
    PURE module function getPrsWNO_XY_RK3(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK2_ENABLED
    PURE module function getPrsWNO_XY_RK2(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK1_ENABLED
    PURE module function getPrsWNO_XY_RK1(x, y) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCor

    ! integer weight: XY.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPrsWTI_XY_CK5(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK4_ENABLED
    PURE module function getPrsWTI_XY_CK4(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK3_ENABLED
    PURE module function getPrsWTI_XY_CK3(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK2_ENABLED
    PURE module function getPrsWTI_XY_CK2(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK1_ENABLED
    PURE module function getPrsWTI_XY_CK1(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPrsWTI_XY_RK5(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK4_ENABLED
    PURE module function getPrsWTI_XY_RK4(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK3_ENABLED
    PURE module function getPrsWTI_XY_RK3(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK2_ENABLED
    PURE module function getPrsWTI_XY_RK2(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK1_ENABLED
    PURE module function getPrsWTI_XY_RK1(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCor

    ! real weight: XY.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPrsWTR_XY_CK5(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK4_ENABLED
    PURE module function getPrsWTR_XY_CK4(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK3_ENABLED
    PURE module function getPrsWTR_XY_CK3(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK2_ENABLED
    PURE module function getPrsWTR_XY_CK2(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

#if CK1_ENABLED
    PURE module function getPrsWTR_XY_CK1(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cor
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPrsWTR_XY_RK5(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK4_ENABLED
    PURE module function getPrsWTR_XY_RK4(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK3_ENABLED
    PURE module function getPrsWTR_XY_RK3(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK2_ENABLED
    PURE module function getPrsWTR_XY_RK2(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

#if RK1_ENABLED
    PURE module function getPrsWTR_XY_RK1(x, y, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cor
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! no weight: ULD.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPrsWNO_ULD_CK5(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                            :: dim
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPrsWNO_ULD_CK4(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                            :: dim
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPrsWNO_ULD_CK3(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                            :: dim
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPrsWNO_ULD_CK2(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                            :: dim
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPrsWNO_ULD_CK1(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                            :: dim
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPrsWNO_ULD_RK5(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPrsWNO_ULD_RK4(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPrsWNO_ULD_RK3(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPrsWNO_ULD_RK2(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPrsWNO_ULD_RK1(sample, dim) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWNO_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCor

    ! integer weight: ULD.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPrsWTI_ULD_CK5(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPrsWTI_ULD_CK4(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPrsWTI_ULD_CK3(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPrsWTI_ULD_CK2(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPrsWTI_ULD_CK1(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPrsWTI_ULD_RK5(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPrsWTI_ULD_RK4(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPrsWTI_ULD_RK3(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPrsWTI_ULD_RK2(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPrsWTI_ULD_RK1(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTI_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                            :: dim
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCor

    ! real weight: ULD.

    interface getCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPrsWTR_ULD_CK5(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPrsWTR_ULD_CK4(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPrsWTR_ULD_CK3(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPrsWTR_ULD_CK2(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPrsWTR_ULD_CK1(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPrsWTR_ULD_RK5(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPrsWTR_ULD_RK4(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPrsWTR_ULD_RK3(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPrsWTR_ULD_RK2(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPrsWTR_ULD_RK1(sample, dim, weight) result(cor)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPrsWTR_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                            :: dim
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the (weighted) correlation matrix corresponding to the input (weighted) covariance matrix or return the (weighted) sample Pearson correlation matrix
    !>  of the input array of shape `(ndim, nsam)` or `(nsam, ndim)` or the Pearson correlation coefficient a pair of (weighted) time series `x(1:nsam)` and `y(1:nsam)`
    !>  where `ndim` is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  This generic interface performs one of the following computational tasks:<br>
    !>  <ol>
    !>      <li>    Compute the correlation matrix corresponding to an input covariance matrix in arbitrary `ndim` dimensions.<br>
    !>      <li>    Compute the Pearson correlation coefficient corresponding to an input pair of time series `x` and `y` of `nsam` observations.<br>
    !>      <li>    Compute the Pearson correlation matrix corresponding to an input multivariate sample of `nsam` observations each with `ndim` attributes.<br>
    !>  </ol>
    !>  See the documentation of the parent module [pm_sampleCor](@ref pm_sampleCor) for algorithmic details and sample correlation matrix definition.<br>
    !>
    !>  \note
    !>  The **sample (Pearson) correlation matrix** can also be readily computed from the sample covariance matrix as in the following example,
    !>  \code{.F90}
    !>
    !>      ! unweighted sample.
    !>
    !>      call setCov(cor, subset, .false., sample, dim)
    !>      call setCor(cor, subset, cor, subset)
    !>
    !>      ! weighted sample.
    !>
    !>      call setCov(cor, subset, .false., sample, dim, weight)
    !>      call setCor(cor, subset, cor, subset)
    !>      !
    !>  \endcode
    !>
    !>  \param[inout]   cor         :   The output or input/output positive semi-definite scalar or square matrix of shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the (Pearson) correlation coefficient or matrix corresponding to the input `sample` or time series `x` and `y` or correlation matrix `cov`, whichever is present.<br>
    !>                                  <ol>
    !>                                      <li>    If `x(:)` and `y(:)` are present, then `cor` shall be a scalar, corresponding to the element \f$r_xy\f$ (`cor(1,2)`) of the computed correlation matrix.<br>
    !>                                              Since the correlation matrix is Hermitian, `cor(2,1) = conjg(cor(1,2)` while the diagonal elements are unit values.<br>
    !>                                      <li>    If `cov` is present, then `cor` shall be a square matrix of shape `shape(cov)`.<br>
    !>                                              On output, only the specified input `subset` will be overwritten with the correlation matrix.<br>
    !>                                              Any elements not in the specified input `subset` remains intact.<br>
    !>                                      <li>    If `sample` is present, then `cor` shall be a square matrix of shape `[size(sample, 3 - dim), size(sample, 3 - dim)]`.<br>
    !>                                              On output, only the specified input `subset` will be overwritten with the correlation matrix.<br>
    !>                                              Any elements not in the specified input `subset` remains intact.<br>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar constant argument that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [low](@ref pm_matrixSubset::low), implying that only the lower subset of the output correlation matrix must be computed.<br>
    !>                                              Specifying this value is useful when the contents of the diagonal elements of the input `cor` must be preserved
    !>                                              because they contain valuable information, e.g., [sample variance](@ref pm_sampleCov).<br>
    !>                                              This option is currently available only if the input argument `cov` is present.<br>
    !>                                              The diagonal elements of `cor` will remain untouched upon return.<br>
    !>                                      <li>    The constant [upp](@ref pm_matrixSubset::upp), implying that only the upper subset of the output correlation matrix must be computed.<br>
    !>                                              Specifying this value is useful when the contents of the diagonal elements of the input `cor` must be preserved
    !>                                              because they contain valuable information, e.g., [sample variance](@ref pm_sampleCov).<br>
    !>                                              This option is currently available only if the input argument `cov` is present.<br>
    !>                                              The diagonal elements of `cor` will remain untouched upon return.<br>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the output correlation matrix must be computed.<br>
    !>                                              This option is available only if either of the input arguments `cov` or `sample` are present.<br>
    !>                                              By definition, all diagonal elements of `cor` will be set to `1`.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the output correlation matrix must be computed.<br>
    !>                                              This option is available only if either of the input arguments `cov` or `sample` are present.<br>
    !>                                              By definition, all diagonal elements of `cor` will be set to `1`.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>  \param[in]      cov         :   The input positive semi-definite square matrix of shape `(1 : ndim, 1 : ndim)` of the same type and kind as the input `cor`,
    !>                                  representing the covariance matrix based upon which the output correlation matrix `cor` is to be computed.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `sample`, `x`, and `y` are missing.)
    !>  \param[in]      subsetv     :   The input scalar constant argument that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input covariance matrix must be used.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input covariance matrix must be used.<br>
    !>                                      <li>    The constant [low](@ref pm_matrixSubset::low), implying that only the lower subset of the input covariance matrix must be used.<br>
    !>                                              This option is available **if and only if** the input argument `stdinv` is present (otherwise, how do we know the variances?).<br>
    !>                                      <li>    The constant [upp](@ref pm_matrixSubset::upp), implying that only the upper subset of the input covariance matrix must be used.<br>
    !>                                              This option is available **if and only if** the input argument `stdinv` is present (otherwise, how do we know the variances?).<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `cov` is present.)
    !>  \param[in]      stdinv      :   The input positive vector of shape `(1 : ndim)` of type `real` of the same kind as the input `cor`,
    !>                                  containing the inverse of the standard deviations of the `ndim` data attributes based upon which the output correlation matrix `cor` is to be computed.<br>
    !>                                  \f{equation}{
    !>                                      \ms{stdinv}(i) = \frac{1} {\sqrt{\ms{cov}(i,i)}} ~,
    !>                                  \f}
    !>                                  (**optional**. It must be present **if and only if** the input argument `subsetv` is set to either [low](@ref pm_matrixSubset::low) or [upp](@ref pm_matrixSubset::upp).)
    !>  \param[in]      mean        :   The input scalar or `contiguous` vector of shape `(1:ndim)` of the same type and kind as the output `cor` containing the `sample` mean.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `sample` is a 1D array, then `mean` must be a scalar.<br>
    !>                                      <li>    If the input `sample` is array, then `mean` must be a vector of size `ndim = size(sample, 3 - dim)`
    !>                                              (i.e., computed along the specified input dimension `dim`).<br>
    !>                                  </ol>
    !>                                  The mean vector can be readily computed via [getMean](@ref pm_sampleMean::getMean) or [setMean](@ref pm_sampleMean::setMean).<br>
    !>                                  (**optional**. default = [getFilled(0., ndim)](@ref pm_arrayFill::getFilled) or `[0., 0.]` depending on the rank of the input `sample` and `dim` or the presence of `x` and `y`.)
    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `cov` and `sample` are missing.)
    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `cov` and `sample` are missing.)
    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `cor`,
    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `cor` must be computed (i.e., the direction of individual observations).<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `cov`, `x`, and `y` are missing.)
    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the correlation matrix must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present and is of rank `2`.)
    !>  \param[in]      weight      :   The `contiguous` vector of length `nsam` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, or
    !>                                      <li>    type `real` of the same kind as the kind of the output `cor`,
    !>                                  </ol>
    !>                                  containing the corresponding weights of individual `nsam` observations in `sample` or the pair of vectors `x` and `y`.<br>
    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled). It can be present **if and only if** the input arguments `sample` or `x` and `y` are present.)
    !>  \param[in]      weisum      :   The input scalar of the same type and kind as the input `weight` containing `sum(weight)`.<br>
    !>                                  This quantity is a byproduct of computing the mean of a sample and is automatically returned by [setMean](@ref pm_sampleMean::setMean).<br>
    !>                                  (**optional**. It **must** be present **if and only if** both `mean` and `weight` arguments are present **and** the input argument `cov` is missing.)
    !>
    !>  \interface{setCor}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCor, only: setCor
    !>
    !>      ! correlation matrix to covariance matrix.
    !>
    !>      call setCor(cor(1:ndim, 1:ndim), subset, cov(1:ndim, 1:ndim), subsetv) ! subsetv = uppDia, lowDia
    !>      call setCor(cor(1:ndim, 1:ndim), subset, cov(1:ndim, 1:ndim), subsetv, stdinv(1:ndim)) ! subsetv = upp, low
    !>
    !>      ! XY time series Pearson correlation coefficient.
    !>
    !>      call setCor(cor, x(1:nsam), y(1:nsam)) ! Pearson correlation coefficient.
    !>      call setCor(cor, x(1:nsam), y(1:nsam), weight(1:nsam)) ! Pearson correlation coefficient.
    !>      call setCor(cor, mean(1:2), x(1:nsam), y(1:nsam)) ! Pearson correlation coefficient.
    !>      call setCor(cor, mean(1:2), x(1:nsam), y(1:nsam), weight(1:nsam), weisum) ! Pearson correlation coefficient.
    !>
    !>      ! sample correlation matrix.
    !>
    !>      call setCor(cor(1:ndim, 1:ndim), subset, sample(:,:), dim)
    !>      call setCor(cor(1:ndim, 1:ndim), subset, sample(:,:), dim, weight(1:nsam))
    !>      call setCor(cor(1:ndim, 1:ndim), subset, mean(1:ndim), sample(:,:), dim)
    !>      call setCor(cor(1:ndim, 1:ndim), subset, mean(1:ndim), sample(:,:), dim, weight(1:nsam), weisum)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `any(x(1) /= x)` must hold for the corresponding input arguments.<br>
    !>  The condition `any(y(1) /= x)` must hold for the corresponding input arguments.<br>
    !>  The input `sample` must contain at least two unique values per sample attribute.<br>
    !>  The condition `0 < sum(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 < size(sample, dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(sample, 3 - dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(cov) == shape(cor))` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(variance)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(mean)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, dim) == size(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(x) == size(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(x) == size(y)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{setCor}
    !>  \include{lineno} example/pm_sampleCor/setCor/main.F90
    !>  \compilef{setCor}
    !>  \output{setCor}
    !>  \include{lineno} example/pm_sampleCor/setCor/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
    !>
    !>  \final{setCor}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
    !>  \AmirShahmoradi, Wednesday 4:13 AM, August 13, 2016, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! RK CFC: upp/low from upp/low.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RK CFC: uppDia/lowDia from upp/low.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_RK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_RK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_RK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_RK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_RK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RK CFC: upp/low from uppDia/lowDia.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! RK CFC: uppDia/lowDia from uppDia/lowDia.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_RK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_RK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_RK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_RK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_RK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: cov(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! RK XY - Prs - WNO

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_RK5(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_RK4(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_RK3(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_RK2(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_RK1(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_RK5(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_RK4(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_RK3(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_RK2(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_RK1(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! RK XY - Prs - WTI

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_RK5(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_RK4(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_RK3(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_RK2(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_RK1(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_RK5(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_RK4(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_RK3(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_RK2(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_RK1(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! RK XY - Prs - WTR

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_RK5(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_RK4(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_RK3(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_RK2(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_RK1(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_RK5(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_RK4(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_RK3(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_RK2(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_RK1(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! RK UXD - Prs - WNO

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_RK5(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_RK4(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_RK3(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_RK2(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_RK1(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_RK5(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_RK4(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_RK3(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_RK2(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_RK1(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! RK UXD - Prs - WTI

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_RK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_RK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_RK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_RK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_RK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_RK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_RK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_RK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_RK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_RK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! RK UXD - Prs - WTI real

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_RK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_RK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_RK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_RK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_RK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_RK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_RK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_RK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_RK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_RK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! RK XLD - Prs - WNO

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_RK5(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_RK4(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_RK3(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_RK2(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_RK1(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_RK5(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_RK4(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_RK3(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_RK2(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_RK1(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! RK XLD - Prs - WTI

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_RK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_RK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_RK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_RK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_RK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_RK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_RK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_RK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_RK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_RK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! RK XLD - Prs - WTI real

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_RK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_RK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_RK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_RK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_RK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_RK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_RK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_RK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_RK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_RK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CK CFC: upp/low from upp/low.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXX_VUXX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(upp_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXX_VXLX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(upp_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLX_VUXX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(low_type)      , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLX_VXLX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(low_type)      , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CK CFC: uppDia/lowDia from upp/low.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXD_VUXX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXD_VXLX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLD_VUXX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(upp_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_CK5(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_CK4(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_CK3(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_CK2(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLD_VXLX_CK1(cor, subset, cov, subsetv, stdinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(low_type)      , intent(in)                            :: subsetv
        real(TKG)           , intent(in)    , contiguous            :: stdinv(:)
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CK CFC: upp/low from uppDia/lowDia.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXX_VUXD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VUXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(upp_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXX_VXLD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXX_VXLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(upp_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLX_VUXD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VUXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(low_type)      , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLX_VXLD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLX_VXLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(low_type)      , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CK CFC: uppDia/lowDia from uppDia/lowDia.

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXD_VUXD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VUXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RUXD_VXLD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RUXD_VXLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLD_VUXD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VUXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_CK5(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_CK4(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_CK3(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_CK2(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCFC_RXLD_VXLD_CK1(cor, subset, cov, subsetv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCFC_RXLD_VXLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetv
        complex(TKG)        , intent(in)    , contiguous            :: cov(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cor(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CK XY - Prs - WNO

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_CK5(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_CK4(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_CK3(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_CK2(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWNO_XY_CK1(cor, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_CK5(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_CK4(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_CK3(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_CK2(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWNO_XY_CK1(cor, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! CK XY - Prs - WTI

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_CK5(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_CK4(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_CK3(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_CK2(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWTI_XY_CK1(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_CK5(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_CK4(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_CK3(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_CK2(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWTI_XY_CK1(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! CK XY - Prs - WTR

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_CK5(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_CK4(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_CK3(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_CK2(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWTR_XY_CK1(cor, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_CK5(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_CK4(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_CK3(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_CK2(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWTR_XY_CK1(cor, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)                       :: cor
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CK UXD - Prs - WNO

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_CK5(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_CK4(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_CK3(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_CK2(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWNO_UXD_CK1(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_CK5(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_CK4(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_CK3(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_CK2(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWNO_UXD_CK1(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! CK UXD - Prs - WTI

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_CK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_CK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_CK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_CK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWTI_UXD_CK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_CK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_CK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_CK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_CK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWTI_UXD_CK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! CK UXD - Prs - WTI real

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_CK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_CK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_CK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_CK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWTR_UXD_CK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_CK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_CK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_CK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_CK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWTR_UXD_CK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CK XLD - Prs - WNO

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_CK5(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_CK4(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_CK3(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_CK2(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWNO_XLD_CK1(cor, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWNO_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_CK5(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_CK4(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_CK3(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_CK2(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWNO_XLD_CK1(cor, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWNO_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! CK XLD - Prs - WTI

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_CK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_CK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_CK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_CK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWTI_XLD_CK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTI_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_CK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_CK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_CK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_CK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWTI_XLD_CK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTI_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    ! CK XLD - Prs - WTI real

    interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_CK5(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_CK4(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_CK3(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_CK2(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsOrgWTR_XLD_CK1(cor, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsOrgWTR_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_CK5(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_CK4(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_CK3(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_CK2(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPrsAvgWTR_XLD_CK1(cor, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPrsAvgWTR_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cor(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Spearman rank correlation matrix of the input (weighted) sample of shape `(ndim, nsam)` or `(nsam, ndim)`
    !>  or the Spearman rank correlation coefficient a pair of (weighted) time series `x(1:nsam)` and `y(1:nsam)` where `ndim`
    !>  is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  This generic interface performs one of the following computational tasks:<br>
    !>  <ol>
    !>      <li>    Compute the Spearman rank correlation coefficient corresponding to an input pair of time series `x` and `y` of `nsam` observations.<br>
    !>      <li>    Compute the Spearman rank correlation matrix corresponding to an input multivariate sample of `nsam` observations each with `ndim` attributes.<br>
    !>  </ol>
    !>  See the documentation of the parent module [pm_sampleCor](@ref pm_sampleCor) for algorithmic details and sample correlation matrix definition.<br>
    !>
    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL,
    !>                                      <li>    type `integer` of kind \IKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
    !>                                  </ol>
    !>                                  or a scalar of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
    !>                                  </ol>
    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `rho` and `sample` are missing.)
    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the input `x`,
    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `rho` and `sample` are missing.)
    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL,
    !>                                      <li>    type `integer` of kind \IKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
    !>                                  </ol>
    !>                                  or a scalar of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
    !>                                  </ol>
    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `rho` must be computed (i.e., the direction of individual observations).<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `rho`, `x`, and `y` are missing.)
    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the correlation matrix must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present and is of rank `2`.)
    !>  \param[in]      weight      :   The input `contiguous` vector of length `nsam`,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK,
    !>                                      <li>    type `real` of default kind \RK,
    !>                                  </ol>
    !>                                  containing the corresponding weights of individual `nsam` observations in `sample` or the pair of vectors `x` and `y`.<br>
    !>                                  Note that this default \RK kind type parameter requirement on input `weight` of type `real` is unlike the other `pm_sample*` modules of the ParaMonte library.<br>
    !>                                  This requirement is enforced by the default kind type parameter of the output of [setRankFractional](@ref pm_arrayRank::setRankFractional).<br>
    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled).)
    !>
    !>  \return
    !>  `rho`                       :   The output positive semi-definite scalar or square matrix of shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `real` of default kind \RK,
    !>                                  </ol>
    !>                                  containing the (Spearman rank) correlation coefficient or full correlation matrix corresponding to the input `sample` or time series `x` and `y`, whichever is present.<br>
    !>                                  <ol>
    !>                                      <li>    If `x(:)` and `y(:)` are present, then `rho` shall be a scalara scalar of value \f$r_{xy}\f$,
    !>                                              \f{equation}{
    !>                                                  \ms{rho} =
    !>                                                  \begin{bmatrix}
    !>                                                      1 && r_{xy} \\
    !>                                                      r_{yx} && 1 ~.
    !>                                                  \end{bmatrix}
    !>                                              \f}
    !>                                      <li>    If `sample` is present, then `rho` shall be a square matrix of shape `[size(sample, 3 - dim), size(sample, 3 - dim)]`.<br>
    !>                                  </ol>
    !>
    !>  \interface{getRho}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCor, only: getRho
    !>
    !>      ! XY time series Spearman rank correlation coefficient.
    !>
    !>      rho = getRho(x(1:nsam), y(1:nsam)                ) ! Spearman rank correlation coefficient.
    !>      rho = getRho(x(1:nsam), y(1:nsam), weight(1:nsam)) ! Spearman rank correlation coefficient.
    !>
    !>      ! sample Spearman rank correlation matrix.
    !>
    !>      rho(1:ndim, 1:ndim) = getRho(sample(:,:), dim)
    !>      rho(1:ndim, 1:ndim) = getRho(sample(:,:), dim, weight(1:nsam))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All conditions that must hold for the generic interface [setRho](@ref pm_sampleCor::setRho) must equally hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{getRho}
    !>  \include{lineno} example/pm_sampleCor/getRho/main.F90
    !>  \compilef{getRho}
    !>  \output{getRho}
    !>  \include{lineno} example/pm_sampleCor/getRho/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
    !>
    !>  \final{getRho}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! XY - Rho - WNO

    interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWNO_XY_D0_SK5(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWNO_XY_D0_SK4(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWNO_XY_D0_SK3(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWNO_XY_D0_SK2(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWNO_XY_D0_SK1(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWNO_XY_D1_SK5(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWNO_XY_D1_SK4(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWNO_XY_D1_SK3(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWNO_XY_D1_SK2(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWNO_XY_D1_SK1(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRhoWNO_XY_D1_IK5(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK4_ENABLED
    PURE module function getRhoWNO_XY_D1_IK4(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK3_ENABLED
    PURE module function getRhoWNO_XY_D1_IK3(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK2_ENABLED
    PURE module function getRhoWNO_XY_D1_IK2(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK1_ENABLED
    PURE module function getRhoWNO_XY_D1_IK1(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRhoWNO_XY_D1_RK5(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK4_ENABLED
    PURE module function getRhoWNO_XY_D1_RK4(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK3_ENABLED
    PURE module function getRhoWNO_XY_D1_RK3(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK2_ENABLED
    PURE module function getRhoWNO_XY_D1_RK2(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK1_ENABLED
    PURE module function getRhoWNO_XY_D1_RK1(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getRhoWNO_XY_D1_PSSK5(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWNO_XY_D1_PSSK4(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWNO_XY_D1_PSSK3(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWNO_XY_D1_PSSK2(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWNO_XY_D1_PSSK1(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getRhoWNO_XY_D1_BSSK(x, y) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_XY_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRho

    ! XY - Rho - WTI

    interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWTI_XY_D0_SK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTI_XY_D0_SK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTI_XY_D0_SK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTI_XY_D0_SK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTI_XY_D0_SK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWTI_XY_D1_SK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTI_XY_D1_SK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTI_XY_D1_SK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTI_XY_D1_SK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTI_XY_D1_SK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRhoWTI_XY_D1_IK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK4_ENABLED
    PURE module function getRhoWTI_XY_D1_IK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK3_ENABLED
    PURE module function getRhoWTI_XY_D1_IK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK2_ENABLED
    PURE module function getRhoWTI_XY_D1_IK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK1_ENABLED
    PURE module function getRhoWTI_XY_D1_IK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRhoWTI_XY_D1_RK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK4_ENABLED
    PURE module function getRhoWTI_XY_D1_RK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK3_ENABLED
    PURE module function getRhoWTI_XY_D1_RK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK2_ENABLED
    PURE module function getRhoWTI_XY_D1_RK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK1_ENABLED
    PURE module function getRhoWTI_XY_D1_RK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getRhoWTI_XY_D1_PSSK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTI_XY_D1_PSSK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTI_XY_D1_PSSK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTI_XY_D1_PSSK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTI_XY_D1_PSSK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getRhoWTI_XY_D1_BSSK(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_XY_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRho

    ! XY - Rho - WTR

    interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWTR_XY_D0_SK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTR_XY_D0_SK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTR_XY_D0_SK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTR_XY_D0_SK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTR_XY_D0_SK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWTR_XY_D1_SK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTR_XY_D1_SK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTR_XY_D1_SK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTR_XY_D1_SK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTR_XY_D1_SK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRhoWTR_XY_D1_IK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK4_ENABLED
    PURE module function getRhoWTR_XY_D1_IK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK3_ENABLED
    PURE module function getRhoWTR_XY_D1_IK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK2_ENABLED
    PURE module function getRhoWTR_XY_D1_IK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if IK1_ENABLED
    PURE module function getRhoWTR_XY_D1_IK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRhoWTR_XY_D1_RK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK4_ENABLED
    PURE module function getRhoWTR_XY_D1_RK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK3_ENABLED
    PURE module function getRhoWTR_XY_D1_RK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK2_ENABLED
    PURE module function getRhoWTR_XY_D1_RK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if RK1_ENABLED
    PURE module function getRhoWTR_XY_D1_RK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getRhoWTR_XY_D1_PSSK5(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTR_XY_D1_PSSK4(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTR_XY_D1_PSSK3(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTR_XY_D1_PSSK2(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTR_XY_D1_PSSK1(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getRhoWTR_XY_D1_BSSK(x, y, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_XY_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)                                                   :: rho
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ULD - Rho - WNO

    interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWNO_ULD_SK5(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWNO_ULD_SK4(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWNO_ULD_SK3(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWNO_ULD_SK2(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWNO_ULD_SK1(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRhoWNO_ULD_IK5(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getRhoWNO_ULD_IK4(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getRhoWNO_ULD_IK3(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getRhoWNO_ULD_IK2(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getRhoWNO_ULD_IK1(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRhoWNO_ULD_RK5(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getRhoWNO_ULD_RK4(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getRhoWNO_ULD_RK3(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getRhoWNO_ULD_RK2(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getRhoWNO_ULD_RK1(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getRhoWNO_ULD_PSSK5(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWNO_ULD_PSSK4(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWNO_ULD_PSSK3(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWNO_ULD_PSSK2(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWNO_ULD_PSSK1(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getRhoWNO_ULD_BSSK(sample, dim) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWNO_ULD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRho

    ! ULD - Rho - WTI

    interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWTI_ULD_SK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTI_ULD_SK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTI_ULD_SK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTI_ULD_SK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTI_ULD_SK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRhoWTI_ULD_IK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getRhoWTI_ULD_IK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getRhoWTI_ULD_IK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getRhoWTI_ULD_IK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getRhoWTI_ULD_IK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRhoWTI_ULD_RK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getRhoWTI_ULD_RK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getRhoWTI_ULD_RK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getRhoWTI_ULD_RK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getRhoWTI_ULD_RK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getRhoWTI_ULD_PSSK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTI_ULD_PSSK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTI_ULD_PSSK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTI_ULD_PSSK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTI_ULD_PSSK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getRhoWTI_ULD_BSSK(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTI_ULD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRho

    ! ULD - Rho - WTR

    interface getRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRhoWTR_ULD_SK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTR_ULD_SK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTR_ULD_SK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTR_ULD_SK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTR_ULD_SK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRhoWTR_ULD_IK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getRhoWTR_ULD_IK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getRhoWTR_ULD_IK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getRhoWTR_ULD_IK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getRhoWTR_ULD_IK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRhoWTR_ULD_RK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getRhoWTR_ULD_RK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getRhoWTR_ULD_RK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getRhoWTR_ULD_RK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getRhoWTR_ULD_RK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getRhoWTR_ULD_PSSK5(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getRhoWTR_ULD_PSSK4(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getRhoWTR_ULD_PSSK3(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getRhoWTR_ULD_PSSK2(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getRhoWTR_ULD_PSSK1(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getRhoWTR_ULD_BSSK(sample, dim, weight) result(rho)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRhoWTR_ULD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        integer(IK)             , intent(in)                        :: dim
        real(TKR)                                                   :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRho

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Spearman rank correlation matrix of the input (weighted) sample of shape `(ndim, nsam)` or `(nsam, ndim)`
    !>  or the Spearman rank correlation coefficient a pair of (weighted) time series `x(1:nsam)` and `y(1:nsam)` where `ndim`
    !>  is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  This generic interface performs one of the following computational tasks:<br>
    !>  <ol>
    !>      <li>    Compute the Spearman rank correlation coefficient corresponding to an input pair of time series `x` and `y` of `nsam` observations.<br>
    !>      <li>    Compute the Spearman rank correlation matrix corresponding to an input multivariate sample of `nsam` observations each with `ndim` attributes.<br>
    !>  </ol>
    !>  See the documentation of the parent module [pm_sampleCor](@ref pm_sampleCor) for algorithmic details and sample correlation matrix definition.<br>
    !>
    !>  \param[inout]   rho         :   The output or input/output positive semi-definite scalar or square matrix of shape `(1 : ndim, 1 : ndim)` of
    !>                                  <ol>
    !>                                      <li>    type `real` of default kind \RK,
    !>                                  </ol>
    !>                                  containing the (Spearman rank) correlation coefficient or matrix corresponding to the input `sample` or time series `x` and `y`, whichever is present.<br>
    !>                                  <ol>
    !>                                      <li>    If `x(:)` and `y(:)` are present, then `rho` shall be a scalar.<br>
    !>                                      <li>    If `sample` is present, then `rho` shall be a square matrix of shape `[size(sample, 3 - dim), size(sample, 3 - dim)]`.<br>
    !>                                              On output, only the specified input `subset` will be overwritten with the correlation matrix.<br>
    !>                                              Any elements not in the specified input `subset` remains intact.<br>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar constant argument that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the output correlation matrix must be computed.<br>
    !>                                              This option is available only if either of the input argument `sample` is present.<br>
    !>                                              By definition, all diagonal elements of `rho` will be set to `1`.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the output correlation matrix must be computed.<br>
    !>                                              This option is available only if either of the input argument `sample` is present.<br>
    !>                                              By definition, all diagonal elements of `rho` will be set to `1`.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>  \param[out]     frankx      :   The output `contiguous` vector of shape `(nsam)` of the same type and kind as the output `rho`,
    !>                                  containing the fractional ranking of the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present)
    !>  \param[out]     franky      :   The output `contiguous` vector of shape `(nsam)` of the same type and kind as the output `rho`,
    !>                                  containing the fractional ranking of the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present.)
    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL,
    !>                                      <li>    type `integer` of kind \IKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
    !>                                  </ol>
    !>                                  or a scalar of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
    !>                                  </ol>
    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `rho` and `sample` are missing.)
    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the input `x`,
    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `rho` and `sample` are missing.)
    !>  \param[out]     frank       :   The output `contiguous` array of the same type and kind as the output `rho` and of the same shape as the input `sample`,
    !>                                  containing the fractional ranking of the `sample` attributes along the specified `dim`.<br>
    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `rho` must be computed (i.e., the direction of individual observations).<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `sample` is present.)
    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL,
    !>                                      <li>    type `integer` of kind \IKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
    !>                                  </ol>
    !>                                  or a scalar of,
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
    !>                                  </ol>
    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `rho` must be computed (i.e., the direction of individual observations).<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `rho`, `x`, and `y` are missing.)
    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the correlation matrix must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present and is of rank `2`.)
    !>  \param[in]      weight      :   The input `contiguous` vector of length `nsam`,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK,
    !>                                      <li>    type `real` of default kind \RK,
    !>                                  </ol>
    !>                                  containing the corresponding weights of individual `nsam` observations in `sample` or the pair of vectors `x` and `y`.<br>
    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled).)
    !>
    !>  \interface{setRho}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCor, only: setRho
    !>
    !>      ! XY time series Spearman rank correlation coefficient.
    !>
    !>      call setRho(rho, frankx(1:nsam), franky(1:nsam), x(1:nsam), y(1:nsam)                ) ! Spearman rank correlation coefficient.
    !>      call setRho(rho, frankx(1:nsam), franky(1:nsam), x(1:nsam), y(1:nsam), weight(1:nsam)) ! Spearman rank correlation coefficient.
    !>
    !>      ! sample Spearman rank correlation matrix.
    !>
    !>      call setRho(rho(1:ndim, 1:ndim), subset, sample(:,:), dim)
    !>      call setRho(rho(1:ndim, 1:ndim), subset, sample(:,:), dim, weight(1:nsam))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `any(x(1) /= x)` must hold for the corresponding input arguments.<br>
    !>  The condition `any(y(1) /= x)` must hold for the corresponding input arguments.<br>
    !>  The input `sample` must contain at least two unique values per sample attribute.<br>
    !>  The condition `0 < sum(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 < size(sample, dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(frankx) == size(x)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(franky) == size(y)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(sample, 3 - dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(frank) == shape(sample))` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(variance)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(mean)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, dim) == size(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(x) == size(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(x) == size(y)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getRho](@ref pm_sampleCor::getRho)<br>
    !>  [setRho](@ref pm_sampleCor::setRho)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{setRho}
    !>  \include{lineno} example/pm_sampleCor/setRho/main.F90
    !>  \compilef{setRho}
    !>  \output{setRho}
    !>  \include{lineno} example/pm_sampleCor/setRho/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
    !>
    !>  \final{setRho}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! XY - Rho - WNO

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_XY_D0_SK5(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_XY_D0_SK4(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_XY_D0_SK3(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_XY_D0_SK2(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_XY_D0_SK1(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_SK5(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_SK4(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_SK3(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_SK2(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_SK1(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_IK5(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_IK4(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_IK3(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_IK2(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_IK1(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_RK5(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_RK5
#endif
        use pm_kind, only: TKG => RK, TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_RK4(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_RK3(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_RK2(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_RK1(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_PSSK5(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_PSSK4(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_PSSK3(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_PSSK2(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_XY_D1_PSSK1(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWNO_XY_D1_BSSK(rho, frankx, franky, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XY_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    ! XY - Rho - WTI

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_XY_D0_SK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_XY_D0_SK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_XY_D0_SK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_XY_D0_SK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_XY_D0_SK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_SK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_SK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_SK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_SK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_SK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_IK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_IK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_IK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_IK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_IK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_RK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_RK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_RK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_RK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_RK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_PSSK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_PSSK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_PSSK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_PSSK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_XY_D1_PSSK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWTI_XY_D1_BSSK(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XY_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    ! XY - Rho - WTR

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_XY_D0_SK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D0_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_XY_D0_SK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D0_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_XY_D0_SK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D0_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_XY_D0_SK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D0_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_XY_D0_SK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D0_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)                        :: x, y
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_SK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_SK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_SK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_SK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_SK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_IK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_IK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_IK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_IK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_IK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_RK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_RK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_RK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_RK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_RK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_PSSK5(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_PSSK4(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_PSSK3(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_PSSK2(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_XY_D1_PSSK1(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWTR_XY_D1_BSSK(rho, frankx, franky, x, y, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XY_D1_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
        real(TKR)               , intent(out)   , contiguous        :: frankx(:), franky(:)
        real(TKR)               , intent(out)                       :: rho
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! UXD - Rho - WNO

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_UXD_SK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_UXD_SK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_UXD_SK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_UXD_SK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_UXD_SK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWNO_UXD_IK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWNO_UXD_IK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWNO_UXD_IK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWNO_UXD_IK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWNO_UXD_IK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWNO_UXD_RK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWNO_UXD_RK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWNO_UXD_RK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWNO_UXD_RK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWNO_UXD_RK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_UXD_PSSK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_UXD_PSSK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_UXD_PSSK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_UXD_PSSK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_UXD_PSSK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWNO_UXD_BSSK(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_UXD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    ! UXD - Rho - WTI

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_UXD_SK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_UXD_SK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_UXD_SK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_UXD_SK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_UXD_SK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWTI_UXD_IK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWTI_UXD_IK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWTI_UXD_IK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWTI_UXD_IK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWTI_UXD_IK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWTI_UXD_RK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWTI_UXD_RK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWTI_UXD_RK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWTI_UXD_RK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWTI_UXD_RK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_UXD_PSSK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_UXD_PSSK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_UXD_PSSK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_UXD_PSSK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_UXD_PSSK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWTI_UXD_BSSK(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_UXD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    ! UXD - Rho - WTR

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_UXD_SK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_UXD_SK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_UXD_SK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_UXD_SK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_UXD_SK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWTR_UXD_IK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWTR_UXD_IK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWTR_UXD_IK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWTR_UXD_IK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWTR_UXD_IK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWTR_UXD_RK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWTR_UXD_RK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWTR_UXD_RK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWTR_UXD_RK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWTR_UXD_RK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_UXD_PSSK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_UXD_PSSK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_UXD_PSSK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_UXD_PSSK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_UXD_PSSK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWTR_UXD_BSSK(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_UXD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! XLD - Rho - WNO

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_XLD_SK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_XLD_SK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_XLD_SK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_XLD_SK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_XLD_SK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWNO_XLD_IK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWNO_XLD_IK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWNO_XLD_IK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWNO_XLD_IK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWNO_XLD_IK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWNO_XLD_RK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWNO_XLD_RK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWNO_XLD_RK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWNO_XLD_RK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWNO_XLD_RK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWNO_XLD_PSSK5(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWNO_XLD_PSSK4(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWNO_XLD_PSSK3(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWNO_XLD_PSSK2(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWNO_XLD_PSSK1(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWNO_XLD_BSSK(rho, subset, frank, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWNO_XLD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    ! XLD - Rho - WTI

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_XLD_SK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_XLD_SK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_XLD_SK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_XLD_SK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_XLD_SK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWTI_XLD_IK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWTI_XLD_IK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWTI_XLD_IK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWTI_XLD_IK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWTI_XLD_IK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWTI_XLD_RK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWTI_XLD_RK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWTI_XLD_RK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWTI_XLD_RK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWTI_XLD_RK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWTI_XLD_PSSK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTI_XLD_PSSK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTI_XLD_PSSK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTI_XLD_PSSK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTI_XLD_PSSK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWTI_XLD_BSSK(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTI_XLD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

    ! XLD - Rho - WTR

    interface setRho

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_XLD_SK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_SK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_XLD_SK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_SK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_XLD_SK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_SK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_XLD_SK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_SK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_XLD_SK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_SK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRhoWTR_XLD_IK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_IK5
#endif
        use pm_kind, only: TKR => RK, IKG => IK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRhoWTR_XLD_IK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_IK4
#endif
        use pm_kind, only: TKR => RK, IKG => IK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRhoWTR_XLD_IK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_IK3
#endif
        use pm_kind, only: TKR => RK, IKG => IK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRhoWTR_XLD_IK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_IK2
#endif
        use pm_kind, only: TKR => RK, IKG => IK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRhoWTR_XLD_IK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_IK1
#endif
        use pm_kind, only: TKR => RK, IKG => IK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRhoWTR_XLD_RK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_RK5
#endif
        use pm_kind, only: TKR => RK, TKG => RK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRhoWTR_XLD_RK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_RK4
#endif
        use pm_kind, only: TKR => RK, TKG => RK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRhoWTR_XLD_RK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_RK3
#endif
        use pm_kind, only: TKR => RK, TKG => RK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRhoWTR_XLD_RK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_RK2
#endif
        use pm_kind, only: TKR => RK, TKG => RK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRhoWTR_XLD_RK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_RK1
#endif
        use pm_kind, only: TKR => RK, TKG => RK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setRhoWTR_XLD_PSSK5(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_PSSK5
#endif
        use pm_kind, only: TKR => RK, SKG => SK5
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRhoWTR_XLD_PSSK4(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_PSSK4
#endif
        use pm_kind, only: TKR => RK, SKG => SK4
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRhoWTR_XLD_PSSK3(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_PSSK3
#endif
        use pm_kind, only: TKR => RK, SKG => SK3
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRhoWTR_XLD_PSSK2(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_PSSK2
#endif
        use pm_kind, only: TKR => RK, SKG => SK2
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRhoWTR_XLD_PSSK1(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_PSSK1
#endif
        use pm_kind, only: TKR => RK, SKG => SK1
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setRhoWTR_XLD_BSSK(rho, subset, frank, sample, dim, weight)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRhoWTR_XLD_BSSK
#endif
        use pm_kind, only: TKR => RK, SKG => SK
        real(TKR)               , intent(in)    , contiguous        :: weight(:)
        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
        real(TKR)               , intent(out)   , contiguous        :: frank(:,:)
        real(TKR)               , intent(inout) , contiguous        :: rho(:,:)
        integer(IK)             , intent(in)                        :: dim
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRho

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Generate and return the Kendall-A rank correlation matrix of the input (weighted) sample of shape `(ndim, nsam)` or `(nsam, ndim)`
!    !>  or the Kendall-A rank correlation coefficient a pair of (weighted) time series `x(1:nsam)` and `y(1:nsam)` where `ndim`
!    !>  is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
!    !>
!    !>  \details
!    !>  This generic interface performs one of the following computational tasks:<br>
!    !>  <ol>
!    !>      <li>    Compute the Kendall-A rank correlation coefficient corresponding to an input pair of time series `x` and `y` of `nsam` observations.<br>
!    !>      <li>    Compute the Kendall-A rank correlation matrix corresponding to an input multivariate sample of `nsam` observations each with `ndim` attributes.<br>
!    !>  </ol>
!    !>  See the documentation of the parent module [pm_sampleCor](@ref pm_sampleCor) for algorithmic details and sample correlation matrix definition.<br>
!    !>
!    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL,
!    !>                                      <li>    type `integer` of kind \IKALL,
!    !>                                      <li>    type `real` of kind \RKALL,
!    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
!    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
!    !>                                  </ol>
!    !>                                  or a scalar of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
!    !>                                  </ol>
!    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `tau` and `sample` are missing.)
!    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the input `x`,
!    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `tau` and `sample` are missing.)
!    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL,
!    !>                                      <li>    type `integer` of kind \IKALL,
!    !>                                      <li>    type `real` of kind \RKALL,
!    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
!    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
!    !>                                  </ol>
!    !>                                  or a scalar of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
!    !>                                  </ol>
!    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
!    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `tau` must be computed (i.e., the direction of individual observations).<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `tau`, `x`, and `y` are missing.)
!    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the correlation matrix must be computed.<br>
!    !>                                  <ol>
!    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
!    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
!    !>                                  </ol>
!    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present and is of rank `2`.)
!    !>  \param[in]      weight      :   The `contiguous` vector of length `nsam` of,
!    !>                                  <ol>
!    !>                                      <li>    type `integer` of default kind \IK, or
!    !>                                      <li>    type `real` of default kind \RK,
!    !>                                  </ol>
!    !>                                  containing the corresponding weights of individual `nsam` observations in `sample` or the pair of vectors `x` and `y`.<br>
!    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled).)
!    !>
!    !>  `tau`                       :   The output positive semi-definite scalar or square matrix of shape `(1 : ndim, 1 : ndim)` of,
!    !>                                  <ol>
!    !>                                      <li>    type `real` of default kind \RK,
!    !>                                  </ol>
!    !>                                  containing the (Kendall-A rank) correlation coefficient or full correlation matrix corresponding to the input `sample` or time series `x` and `y`, whichever is present.<br>
!    !>                                  <ol>
!    !>                                      <li>    If `x(:)` and `y(:)` are present, then `tau` shall be a scalar.<br>
!    !>                                      <li>    If `sample` is present, then `tau` shall be a square matrix of shape `[size(sample, 3 - dim), size(sample, 3 - dim)]`.<br>
!    !>                                  </ol>
!    !>
!    !>  \interface{getTau}
!    !>  \code{.F90}
!    !>
!    !>      use pm_sampleCor, only: getTau
!    !>
!    !>      ! XY time series Kendall-A rank correlation coefficient.
!    !>
!    !>      tau = getTau(x(1:nsam), y(1:nsam)                ) ! Kendall-A rank correlation coefficient.
!    !>      tau = getTau(x(1:nsam), y(1:nsam), weight(1:nsam)) ! Kendall-A rank correlation coefficient.
!    !>
!    !>      ! sample Kendall-A rank correlation matrix.
!    !>
!    !>      tau(1:ndim, 1:ndim) = getTau(sample(:,:), dim)
!    !>      tau(1:ndim, 1:ndim) = getTau(sample(:,:), dim, weight(1:nsam))
!    !>
!    !>  \endcode
!    !>
!    !>  \warning
!    !>  All conditions that must hold for the generic interface [setTau](@ref pm_sampleCor::setTau) must equally hold for this generic interface.<br>
!    !>  \vericons
!    !>
!    !>  \warnpure
!    !>
!    !>  \see
!    !>  [getCor](@ref pm_sampleCor::getCor)<br>
!    !>  [setCor](@ref pm_sampleCor::setCor)<br>
!    !>  [getRho](@ref pm_sampleCor::getRho)<br>
!    !>  [setRho](@ref pm_sampleCor::setRho)<br>
!    !>  [getTau](@ref pm_sampleCor::getTau)<br>
!    !>  [setTau](@ref pm_sampleCor::setTau)<br>
!    !>  [getCov](@ref pm_sampleCov::getCov)<br>
!    !>  [setCov](@ref pm_sampleCov::setCov)<br>
!    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
!    !>  [getMean](@ref pm_sampleMean::getMean)<br>
!    !>  [setMean](@ref pm_sampleMean::setMean)<br>
!    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
!    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
!    !>  [getVar](@ref pm_sampleVar::getVar)<br>
!    !>  [setVar](@ref pm_sampleVar::setVar)<br>
!    !>
!    !>  \example{getTau}
!    !>  \include{lineno} example/pm_sampleCor/getTau/main.F90
!    !>  \compilef{getTau}
!    !>  \output{getTau}
!    !>  \include{lineno} example/pm_sampleCor/getTau/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
!    !>
!    !>  \final{getTau}
!    !>
!    !>  \author
!    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
!    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
!
!    ! XY - Tau - WNO
!
!    interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWNO_XY_D0_SK5(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D0_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWNO_XY_D0_SK4(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D0_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWNO_XY_D0_SK3(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D0_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWNO_XY_D0_SK2(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D0_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWNO_XY_D0_SK1(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D0_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWNO_XY_D1_SK5(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWNO_XY_D1_SK4(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWNO_XY_D1_SK3(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWNO_XY_D1_SK2(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWNO_XY_D1_SK1(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getTauWNO_XY_D1_IK5(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getTauWNO_XY_D1_IK4(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getTauWNO_XY_D1_IK3(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getTauWNO_XY_D1_IK2(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getTauWNO_XY_D1_IK1(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module function getTauWNO_XY_D1_RK5(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE module function getTauWNO_XY_D1_RK4(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE module function getTauWNO_XY_D1_RK3(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE module function getTauWNO_XY_D1_RK2(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE module function getTauWNO_XY_D1_RK1(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module function getTauWNO_XY_D1_PSSK5(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWNO_XY_D1_PSSK4(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWNO_XY_D1_PSSK3(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWNO_XY_D1_PSSK2(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWNO_XY_D1_PSSK1(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module function getTauWNO_XY_D1_BSSK(x, y) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_XY_D1_BSSK
!#endif
!        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface getTau
!
!    ! XY - Tau - WTI
!
!    interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWTI_XY_D0_SK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D0_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTI_XY_D0_SK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D0_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTI_XY_D0_SK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D0_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTI_XY_D0_SK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D0_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTI_XY_D0_SK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D0_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWTI_XY_D1_SK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTI_XY_D1_SK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTI_XY_D1_SK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTI_XY_D1_SK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTI_XY_D1_SK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getTauWTI_XY_D1_IK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getTauWTI_XY_D1_IK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getTauWTI_XY_D1_IK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getTauWTI_XY_D1_IK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getTauWTI_XY_D1_IK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module function getTauWTI_XY_D1_RK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE module function getTauWTI_XY_D1_RK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE module function getTauWTI_XY_D1_RK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE module function getTauWTI_XY_D1_RK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE module function getTauWTI_XY_D1_RK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module function getTauWTI_XY_D1_PSSK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTI_XY_D1_PSSK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTI_XY_D1_PSSK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTI_XY_D1_PSSK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTI_XY_D1_PSSK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module function getTauWTI_XY_D1_BSSK(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_XY_D1_BSSK
!#endif
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface getTau
!
!    ! XY - Tau - WTR
!
!    interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWTR_XY_D0_SK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D0_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTR_XY_D0_SK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D0_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTR_XY_D0_SK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D0_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTR_XY_D0_SK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D0_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTR_XY_D0_SK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D0_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWTR_XY_D1_SK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTR_XY_D1_SK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTR_XY_D1_SK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTR_XY_D1_SK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTR_XY_D1_SK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getTauWTR_XY_D1_IK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getTauWTR_XY_D1_IK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getTauWTR_XY_D1_IK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getTauWTR_XY_D1_IK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getTauWTR_XY_D1_IK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module function getTauWTR_XY_D1_RK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE module function getTauWTR_XY_D1_RK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE module function getTauWTR_XY_D1_RK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE module function getTauWTR_XY_D1_RK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE module function getTauWTR_XY_D1_RK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module function getTauWTR_XY_D1_PSSK5(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTR_XY_D1_PSSK4(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTR_XY_D1_PSSK3(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTR_XY_D1_PSSK2(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTR_XY_D1_PSSK1(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module function getTauWTR_XY_D1_BSSK(x, y, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_XY_D1_BSSK
!#endif
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)                                                   :: tau
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    ! ULD - Tau - WNO
!
!    interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWNO_ULD_SK5(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWNO_ULD_SK4(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWNO_ULD_SK3(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWNO_ULD_SK2(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWNO_ULD_SK1(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getTauWNO_ULD_IK5(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getTauWNO_ULD_IK4(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getTauWNO_ULD_IK3(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getTauWNO_ULD_IK2(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getTauWNO_ULD_IK1(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module function getTauWNO_ULD_RK5(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE module function getTauWNO_ULD_RK4(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE module function getTauWNO_ULD_RK3(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE module function getTauWNO_ULD_RK2(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE module function getTauWNO_ULD_RK1(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module function getTauWNO_ULD_PSSK5(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWNO_ULD_PSSK4(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWNO_ULD_PSSK3(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWNO_ULD_PSSK2(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWNO_ULD_PSSK1(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module function getTauWNO_ULD_BSSK(sample, dim) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWNO_ULD_BSSK
!#endif
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface getTau
!
!    ! ULD - Tau - WTI
!
!    interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWTI_ULD_SK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTI_ULD_SK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTI_ULD_SK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTI_ULD_SK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTI_ULD_SK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getTauWTI_ULD_IK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getTauWTI_ULD_IK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getTauWTI_ULD_IK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getTauWTI_ULD_IK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getTauWTI_ULD_IK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module function getTauWTI_ULD_RK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE module function getTauWTI_ULD_RK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE module function getTauWTI_ULD_RK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE module function getTauWTI_ULD_RK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE module function getTauWTI_ULD_RK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module function getTauWTI_ULD_PSSK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTI_ULD_PSSK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTI_ULD_PSSK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTI_ULD_PSSK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTI_ULD_PSSK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module function getTauWTI_ULD_BSSK(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTI_ULD_BSSK
!#endif
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface getTau
!
!    ! ULD - Tau - WTR
!
!    interface getTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module function getTauWTR_ULD_SK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTR_ULD_SK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTR_ULD_SK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTR_ULD_SK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTR_ULD_SK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module function getTauWTR_ULD_IK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK4_ENABLED
!    PURE module function getTauWTR_ULD_IK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK3_ENABLED
!    PURE module function getTauWTR_ULD_IK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK2_ENABLED
!    PURE module function getTauWTR_ULD_IK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if IK1_ENABLED
!    PURE module function getTauWTR_ULD_IK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module function getTauWTR_ULD_RK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE module function getTauWTR_ULD_RK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE module function getTauWTR_ULD_RK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE module function getTauWTR_ULD_RK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE module function getTauWTR_ULD_RK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module function getTauWTR_ULD_PSSK5(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK4_ENABLED
!    PURE module function getTauWTR_ULD_PSSK4(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK3_ENABLED
!    PURE module function getTauWTR_ULD_PSSK3(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK2_ENABLED
!    PURE module function getTauWTR_ULD_PSSK2(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#if SK1_ENABLED
!    PURE module function getTauWTR_ULD_PSSK1(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module function getTauWTR_ULD_BSSK(sample, dim, weight) result(tau)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getTauWTR_ULD_BSSK
!#endif
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        real(TKG)                                                   :: tau(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
!    end function
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface getTau
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !>  Return the Kendall-A rank correlation matrix of the input (weighted) sample of shape `(ndim, nsam)` or `(nsam, ndim)`
!    !>  or the Kendall-A rank correlation coefficient a pair of (weighted) time series `x(1:nsam)` and `y(1:nsam)` where `ndim`
!    !>  is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
!    !>
!    !>  \details
!    !>  This generic interface performs one of the following computational tasks:<br>
!    !>  <ol>
!    !>      <li>    Compute the Kendall-A rank correlation coefficient corresponding to an input pair of time series `x` and `y` of `nsam` observations.<br>
!    !>      <li>    Compute the Kendall-A rank correlation matrix corresponding to an input multivariate sample of `nsam` observations each with `ndim` attributes.<br>
!    !>  </ol>
!    !>  See the documentation of the parent module [pm_sampleCor](@ref pm_sampleCor) for algorithmic details and sample correlation matrix definition.<br>
!    !>
!    !>  \param[inout]   tau         :   The output or input/output positive semi-definite scalar or square matrix of shape `(1 : ndim, 1 : ndim)` of type `real` of default kind \RK,
!    !>                                  containing the (Kendall-A rank) correlation coefficient or matrix corresponding to the input `sample` or time series `x` and `y`, whichever is present.<br>
!    !>                                  <ol>
!    !>                                      <li>    If `x(:)` and `y(:)` are present, then `tau` shall be a scalar.<br>
!    !>                                      <li>    If `sample` is present, then `tau` shall be a square matrix of shape `[size(sample, 3 - dim), size(sample, 3 - dim)]`.<br>
!    !>                                              On output, only the specified input `subset` will be overwritten with the correlation matrix.<br>
!    !>                                              Any elements not in the specified input `subset` remains intact.<br>
!    !>                                  </ol>
!    !>  \param[in]      subset      :   The input scalar constant argument that can be any of the following:<br>
!    !>                                  <ol>
!    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the output correlation matrix must be computed.<br>
!    !>                                              This option is available only if either of the input argument `sample` is present.<br>
!    !>                                              By definition, all diagonal elements of `tau` will be set to `1`.<br>
!    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the output correlation matrix must be computed.<br>
!    !>                                              This option is available only if either of the input argument `sample` is present.<br>
!    !>                                              By definition, all diagonal elements of `tau` will be set to `1`.<br>
!    !>                                  </ol>
!    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
!    !>  \param[out]     frankx      :   The output `contiguous` vector of shape `(nsam)` of the same type and kind as the output `tau`,
!    !>                                  containing the fractional ranking of the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present)
!    !>  \param[out]     franky      :   The output `contiguous` vector of shape `(nsam)` of the same type and kind as the output `tau`,
!    !>                                  containing the fractional ranking of the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present.)
!    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL,
!    !>                                      <li>    type `integer` of kind \IKALL,
!    !>                                      <li>    type `real` of kind \RKALL,
!    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
!    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
!    !>                                  </ol>
!    !>                                  or a scalar of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
!    !>                                  </ol>
!    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `tau` and `sample` are missing.)
!    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the input `x`,
!    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `tau` and `sample` are missing.)
!    !>  \param[out]     frank       :   The output `contiguous` array of the same type and kind as the output `tau` and of the same shape as the input `sample`,
!    !>                                  containing the fractional ranking of the `sample` attributes along the specified `dim`.<br>
!    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `tau` must be computed (i.e., the direction of individual observations).<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input argument `sample` is present.)
!    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL,
!    !>                                      <li>    type `integer` of kind \IKALL,
!    !>                                      <li>    type `real` of kind \RKALL,
!    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
!    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
!    !>                                  </ol>
!    !>                                  or a scalar of,
!    !>                                  <ol>
!    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter,
!    !>                                  </ol>
!    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
!    !>                                  If `sample` is a matrix, then the input argument `dim` dictates the direction along which the correlation matrix `tau` must be computed (i.e., the direction of individual observations).<br>
!    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `tau`, `x`, and `y` are missing.)
!    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the correlation matrix must be computed.<br>
!    !>                                  <ol>
!    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
!    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
!    !>                                  </ol>
!    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present and is of rank `2`.)
!    !>  \param[in]      weight      :   The input `contiguous` vector of length `nsam` of,
!    !>                                  <ol>
!    !>                                      <li>    type `integer` of default kind \IK, or
!    !>                                      <li>    type `real` of default kind \RK,
!    !>                                  </ol>
!    !>                                  containing the corresponding weights of individual `nsam` observations in `sample` or the pair of vectors `x` and `y`.<br>
!    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled).)
!    !>
!    !>  \interface{setTau}
!    !>  \code{.F90}
!    !>
!    !>      use pm_sampleCor, only: setTau
!    !>
!    !>      ! XY time series Kendall-A rank correlation coefficient.
!    !>
!    !>      call setTau(tau, frankx(1:nsam), franky(1:nsam), x(1:nsam), y(1:nsam)                ) ! Kendall-A rank correlation coefficient.
!    !>      call setTau(tau, frankx(1:nsam), franky(1:nsam), x(1:nsam), y(1:nsam), weight(1:nsam)) ! Kendall-A rank correlation coefficient.
!    !>
!    !>      ! sample Kendall-A rank correlation matrix.
!    !>
!    !>      call setTau(tau(1:ndim, 1:ndim), subset, sample(:,:), dim)
!    !>      call setTau(tau(1:ndim, 1:ndim), subset, sample(:,:), dim, weight(1:nsam))
!    !>
!    !>  \endcode
!    !>
!    !>  \warning
!    !>  The condition `0 < sum(weight)` must hold for the corresponding input arguments.<br>
!    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
!    !>  The condition `1 < size(sample, dim)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(frankx) == size(x)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(franky) == size(y)` must hold for the corresponding input arguments.<br>
!    !>  The condition `0 < size(sample, 3 - dim)` must hold for the corresponding input arguments.<br>
!    !>  The condition `all(shape(frank) == shape(sample))` must hold for the corresponding input arguments.<br>
!    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(sample, 3 - dim) == size(variance)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(sample, 3 - dim) == size(mean)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(sample, dim) == size(weight)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(x) == size(weight)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(x) == size(y)` must hold for the corresponding input arguments.<br>
!    !>  \vericons
!    !>
!    !>  \warnpure
!    !>
!    !>  \see
!    !>  [getCor](@ref pm_sampleCor::getCor)<br>
!    !>  [setCor](@ref pm_sampleCor::setCor)<br>
!    !>  [getRho](@ref pm_sampleCor::getRho)<br>
!    !>  [setRho](@ref pm_sampleCor::setRho)<br>
!    !>  [getTau](@ref pm_sampleCor::getTau)<br>
!    !>  [setTau](@ref pm_sampleCor::setTau)<br>
!    !>  [getCov](@ref pm_sampleCov::getCov)<br>
!    !>  [setCov](@ref pm_sampleCov::setCov)<br>
!    !>  [setECDF](@ref pm_sampleECDF::setECDF)<br>
!    !>  [getMean](@ref pm_sampleMean::getMean)<br>
!    !>  [setMean](@ref pm_sampleMean::setMean)<br>
!    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
!    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
!    !>  [getVar](@ref pm_sampleVar::getVar)<br>
!    !>  [setVar](@ref pm_sampleVar::setVar)<br>
!    !>
!    !>  \example{setTau}
!    !>  \include{lineno} example/pm_sampleCor/setTau/main.F90
!    !>  \compilef{setTau}
!    !>  \output{setTau}
!    !>  \include{lineno} example/pm_sampleCor/setTau/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
!    !>
!    !>  \final{setTau}
!    !>
!    !>  \author
!    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>
!    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
!
!    ! XY - Tau - WNO
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_XY_D0_SK5(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D0_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_XY_D0_SK4(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D0_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_XY_D0_SK3(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D0_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_XY_D0_SK2(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D0_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_XY_D0_SK1(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D0_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_SK5(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_SK4(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_SK3(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_SK2(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_SK1(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_IK5(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_IK4(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_IK3(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_IK2(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_IK1(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_RK5(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_RK4(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_RK3(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_RK2(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_RK1(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_PSSK5(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_PSSK4(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_PSSK3(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_PSSK2(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_XY_D1_PSSK1(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWNO_XY_D1_BSSK(tau, frankx, franky, x, y)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XY_D1_BSSK
!#endif
!        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    ! XY - Tau - WTI
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_XY_D0_SK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D0_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_XY_D0_SK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D0_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_XY_D0_SK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D0_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_XY_D0_SK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D0_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_XY_D0_SK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D0_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_SK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_SK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_SK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_SK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_SK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_IK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_IK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_IK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_IK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_IK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_RK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_RK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_RK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_RK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_RK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_PSSK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_PSSK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_PSSK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_PSSK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_XY_D1_PSSK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWTI_XY_D1_BSSK(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XY_D1_BSSK
!#endif
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    ! XY - Tau - WTR
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_XY_D0_SK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D0_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_XY_D0_SK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D0_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_XY_D0_SK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D0_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_XY_D0_SK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D0_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_XY_D0_SK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D0_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)                        :: x, y
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_SK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_SK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_SK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_SK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_SK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_IK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_IK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_IK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_IK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_IK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_RK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_RK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_RK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_RK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_RK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_PSSK5(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_PSSK4(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_PSSK3(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_PSSK2(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_XY_D1_PSSK1(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWTR_XY_D1_BSSK(tau, frankx, franky, x, y, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XY_D1_BSSK
!#endif
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: x(:), y(:)
!        real(TKG)               , intent(out)   , contiguous        :: frankx(:), franky(:)
!        real(TKG)               , intent(out)                       :: tau
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    ! UXD - Tau - WNO
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_UXD_SK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_UXD_SK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_UXD_SK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_UXD_SK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_UXD_SK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWNO_UXD_IK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWNO_UXD_IK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWNO_UXD_IK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWNO_UXD_IK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWNO_UXD_IK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWNO_UXD_RK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWNO_UXD_RK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWNO_UXD_RK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWNO_UXD_RK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWNO_UXD_RK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_UXD_PSSK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_UXD_PSSK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_UXD_PSSK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_UXD_PSSK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_UXD_PSSK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWNO_UXD_BSSK(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_UXD_BSSK
!#endif
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    ! UXD - Tau - WTI
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_UXD_SK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_UXD_SK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_UXD_SK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_UXD_SK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_UXD_SK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWTI_UXD_IK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWTI_UXD_IK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWTI_UXD_IK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWTI_UXD_IK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWTI_UXD_IK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWTI_UXD_RK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWTI_UXD_RK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWTI_UXD_RK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWTI_UXD_RK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWTI_UXD_RK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_UXD_PSSK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_UXD_PSSK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_UXD_PSSK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_UXD_PSSK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_UXD_PSSK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWTI_UXD_BSSK(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_UXD_BSSK
!#endif
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    ! UXD - Tau - WTR
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_UXD_SK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_UXD_SK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_UXD_SK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_UXD_SK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_UXD_SK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWTR_UXD_IK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWTR_UXD_IK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWTR_UXD_IK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWTR_UXD_IK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWTR_UXD_IK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWTR_UXD_RK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWTR_UXD_RK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWTR_UXD_RK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWTR_UXD_RK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWTR_UXD_RK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_UXD_PSSK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_UXD_PSSK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_UXD_PSSK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_UXD_PSSK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_UXD_PSSK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWTR_UXD_BSSK(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_UXD_BSSK
!#endif
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(uppDia_type)       , intent(in)                        :: subset
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    ! XLD - Tau - WNO
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_XLD_SK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_XLD_SK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_XLD_SK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_XLD_SK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_XLD_SK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWNO_XLD_IK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWNO_XLD_IK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWNO_XLD_IK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWNO_XLD_IK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWNO_XLD_IK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWNO_XLD_RK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWNO_XLD_RK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWNO_XLD_RK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWNO_XLD_RK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWNO_XLD_RK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWNO_XLD_PSSK5(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWNO_XLD_PSSK4(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWNO_XLD_PSSK3(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWNO_XLD_PSSK2(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWNO_XLD_PSSK1(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWNO_XLD_BSSK(tau, subset, frank, sample, dim)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWNO_XLD_BSSK
!#endif
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    ! XLD - Tau - WTI
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_XLD_SK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_XLD_SK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_XLD_SK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_XLD_SK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_XLD_SK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWTI_XLD_IK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWTI_XLD_IK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWTI_XLD_IK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWTI_XLD_IK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWTI_XLD_IK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWTI_XLD_RK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWTI_XLD_RK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWTI_XLD_RK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWTI_XLD_RK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWTI_XLD_RK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTI_XLD_PSSK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTI_XLD_PSSK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTI_XLD_PSSK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTI_XLD_PSSK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTI_XLD_PSSK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWTI_XLD_BSSK(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTI_XLD_BSSK
!#endif
!        integer(IK)             , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau
!
!    ! XLD - Tau - WTR
!
!    interface setTau
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_XLD_SK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_SK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_XLD_SK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_SK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_XLD_SK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_SK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_XLD_SK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_SK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_XLD_SK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_SK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        character(*,SKG)        , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setTauWTR_XLD_IK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_IK5
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setTauWTR_XLD_IK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_IK4
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setTauWTR_XLD_IK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_IK3
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setTauWTR_XLD_IK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_IK2
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setTauWTR_XLD_IK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_IK1
!#endif
!        use pm_kind, only: TKG => RK, IKG => IK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        integer(IKG)            , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setTauWTR_XLD_RK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_RK5
!#endif
!        use pm_kind, only: TKG => RK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setTauWTR_XLD_RK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_RK4
!#endif
!        use pm_kind, only: TKG => RK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setTauWTR_XLD_RK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_RK3
!#endif
!        use pm_kind, only: TKG => RK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setTauWTR_XLD_RK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_RK2
!#endif
!        use pm_kind, only: TKG => RK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setTauWTR_XLD_RK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_RK1
!#endif
!        use pm_kind, only: TKG => RK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if PDT_ENABLED
!
!#if SK5_ENABLED
!    PURE module subroutine setTauWTR_XLD_PSSK5(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_PSSK5
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK5
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setTauWTR_XLD_PSSK4(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_PSSK4
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK4
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setTauWTR_XLD_PSSK3(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_PSSK3
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK3
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setTauWTR_XLD_PSSK2(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_PSSK2
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK2
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setTauWTR_XLD_PSSK1(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_PSSK1
!#endif
!        use pm_kind, only: TKG => RK, SKG => SK1
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_pdt(SKG))      , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!#endif
!
!#endif
!PDT_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setTauWTR_XLD_BSSK(tau, subset, frank, sample, dim, weight)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setTauWTR_XLD_BSSK
!#endif
!        real(TKG)               , intent(in)    , contiguous        :: weight(:)
!        type(css_type)          , intent(in)    , contiguous        :: sample(:,:)
!        real(TKG)               , intent(out)   , contiguous        :: frank(:,:)
!        real(TKG)               , intent(inout) , contiguous        :: tau(:,:)
!        integer(IK)             , intent(in)                        :: dim
!        type(lowDia_type)       , intent(in)                        :: subset
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setTau

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return the Cordance vector of the input data series `x` and `y`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleCor](@ref pm_sampleCor) for the definition of sample cordance.<br>
    !>
    !>  \param[out]     concordance :   The output scalar of type `integer` of default kind \IK, containing the number of concordant pairs in the sample.<br>
    !>                                  (**optional**. It must be present **if and only if** the output argument `cordance` is missing.)
    !>  \param[out]     discordance :   The output scalar of type `integer` of default kind \IK, containing the number of discordant pairs in the sample.<br>
    !>                                  (**optional**. It must be present **if and only if** the output argument `cordance` is missing.)
    !>  \param[out]     cordance    :   The output scalar of type `integer` of default kind \IK, containing the number of concordant minus discordant pairs in the sample.<br>
    !>                                  (**optional**. It must be present **if and only if** the output arguments `concordance` and `discordance` are both missing.)
    !>  \param[out]     tiedx       :   The output scalar of type `integer` of default kind \IK, containing the number of pairs whose `x` corresponding values are equal (tied).<br>
    !>  \param[out]     tiedy       :   The output scalar of type `integer` of default kind \IK, containing the number of pairs whose `y` corresponding values are equal (tied).<br>
    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of,<br>
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL,<br>
    !>                                      <li>    type `integer` of kind \IKALL,<br>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                      <li>    type string container [css_type](@ref pm_container::css_type),
    !>                                      <li>    type string PDT container [css_pdt](@ref pm_container::css_pdt),
    !>                                  </ol>
    !>                                  or scalar `character` of kind \SKALL of arbitrary length type parameter,
    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as `x`,
    !>                                  containing the second attribute `y` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>
    !>  \interface{setCordance}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCor, only: setCordance
    !>
    !>      call setCordance(cordance, tiedx, tiedy, x, y)
    !>      call setCordance(concordance, discordance, tiedx, tiedy, x, y)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_sampleCor/setCordance/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_sampleCor/setCordance/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCor](@ref test_pm_sampleCor)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface should be extended to allow custom user-specified cordance criteria.<br>
    !>
    !>  \final{setCordance}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! WNO

    interface setCordance

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCordanceSum_D0_SK5(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D0_SK5
#endif
        use pm_kind, only: TKG => RK, SKG => SK5
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCordanceSum_D0_SK4(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D0_SK4
#endif
        use pm_kind, only: TKG => RK, SKG => SK4
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCordanceSum_D0_SK3(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D0_SK3
#endif
        use pm_kind, only: TKG => RK, SKG => SK3
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCordanceSum_D0_SK2(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D0_SK2
#endif
        use pm_kind, only: TKG => RK, SKG => SK2
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCordanceSum_D0_SK1(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D0_SK1
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCordanceSum_D1_SK5(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_SK5
#endif
        use pm_kind, only: TKG => RK, SKG => SK5
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCordanceSum_D1_SK4(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_SK4
#endif
        use pm_kind, only: TKG => RK, SKG => SK4
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCordanceSum_D1_SK3(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_SK3
#endif
        use pm_kind, only: TKG => RK, SKG => SK3
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCordanceSum_D1_SK2(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_SK2
#endif
        use pm_kind, only: TKG => RK, SKG => SK2
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCordanceSum_D1_SK1(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_SK1
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCordanceSum_D1_IK5(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_IK5
#endif
        use pm_kind, only: TKG => RK, IKG => IK5
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCordanceSum_D1_IK4(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_IK4
#endif
        use pm_kind, only: TKG => RK, IKG => IK4
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCordanceSum_D1_IK3(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_IK3
#endif
        use pm_kind, only: TKG => RK, IKG => IK3
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCordanceSum_D1_IK2(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_IK2
#endif
        use pm_kind, only: TKG => RK, IKG => IK2
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCordanceSum_D1_IK1(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_IK1
#endif
        use pm_kind, only: TKG => RK, IKG => IK1
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCordanceSum_D1_RK5(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCordanceSum_D1_RK4(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCordanceSum_D1_RK3(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCordanceSum_D1_RK2(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCordanceSum_D1_RK1(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setCordanceSum_D1_PSSK5(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_PSSK5
#endif
        use pm_kind, only: TKG => RK, SKG => SK5
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCordanceSum_D1_PSSK4(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_PSSK4
#endif
        use pm_kind, only: TKG => RK, SKG => SK4
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCordanceSum_D1_PSSK3(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_PSSK3
#endif
        use pm_kind, only: TKG => RK, SKG => SK3
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCordanceSum_D1_PSSK2(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_PSSK2
#endif
        use pm_kind, only: TKG => RK, SKG => SK2
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCordanceSum_D1_PSSK1(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_PSSK1
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setCordanceSum_D1_BSSK(cordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceSum_D1_BSSK
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        type(css_type)      , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: cordance
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCordanceAll_D0_SK5(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D0_SK5
#endif
        use pm_kind, only: TKG => RK, SKG => SK5
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCordanceAll_D0_SK4(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D0_SK4
#endif
        use pm_kind, only: TKG => RK, SKG => SK4
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCordanceAll_D0_SK3(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D0_SK3
#endif
        use pm_kind, only: TKG => RK, SKG => SK3
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCordanceAll_D0_SK2(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D0_SK2
#endif
        use pm_kind, only: TKG => RK, SKG => SK2
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCordanceAll_D0_SK1(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D0_SK1
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        character(*,SKG)    , intent(in)                            :: x, y
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setCordanceAll_D1_SK5(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_SK5
#endif
        use pm_kind, only: TKG => RK, SKG => SK5
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCordanceAll_D1_SK4(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_SK4
#endif
        use pm_kind, only: TKG => RK, SKG => SK4
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCordanceAll_D1_SK3(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_SK3
#endif
        use pm_kind, only: TKG => RK, SKG => SK3
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCordanceAll_D1_SK2(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_SK2
#endif
        use pm_kind, only: TKG => RK, SKG => SK2
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCordanceAll_D1_SK1(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_SK1
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        character(*,SKG)    , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setCordanceAll_D1_IK5(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_IK5
#endif
        use pm_kind, only: TKG => RK, IKG => IK5
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setCordanceAll_D1_IK4(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_IK4
#endif
        use pm_kind, only: TKG => RK, IKG => IK4
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setCordanceAll_D1_IK3(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_IK3
#endif
        use pm_kind, only: TKG => RK, IKG => IK3
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setCordanceAll_D1_IK2(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_IK2
#endif
        use pm_kind, only: TKG => RK, IKG => IK2
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setCordanceAll_D1_IK1(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_IK1
#endif
        use pm_kind, only: TKG => RK, IKG => IK1
        integer(IKG)        , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCordanceAll_D1_RK5(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCordanceAll_D1_RK4(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCordanceAll_D1_RK3(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCordanceAll_D1_RK2(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCordanceAll_D1_RK1(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setCordanceAll_D1_PSSK5(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_PSSK5
#endif
        use pm_kind, only: TKG => RK, SKG => SK5
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setCordanceAll_D1_PSSK4(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_PSSK4
#endif
        use pm_kind, only: TKG => RK, SKG => SK4
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setCordanceAll_D1_PSSK3(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_PSSK3
#endif
        use pm_kind, only: TKG => RK, SKG => SK3
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setCordanceAll_D1_PSSK2(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_PSSK2
#endif
        use pm_kind, only: TKG => RK, SKG => SK2
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setCordanceAll_D1_PSSK1(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_PSSK1
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        type(css_pdt(SKG))  , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setCordanceAll_D1_BSSK(concordance, discordance, tiedx, tiedy, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCordanceAll_D1_BSSK
#endif
        use pm_kind, only: TKG => RK, SKG => SK1
        type(css_type)      , intent(in)    , contiguous            :: x(:), y(:)
        integer(IK)         , intent(out)                           :: tiedx, tiedy
        integer(IK)         , intent(out)                           :: concordance, discordance
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleCor ! LCOV_EXCL_LINE
