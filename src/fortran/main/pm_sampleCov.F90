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
!>  This module contains classes and procedures for computing the properties related to the covariance matrices of a random sample.
!>
!>  \details
!>
!>  Covariance
!>  ==========
!>
!>  The concept of [variance](@ref pm_sampleVar) can be generalized to measure the covariation of any pair of data attributes.<br>
!>  The **sample covariance matrix** is a \f$K\f$-by-\f$K\f$ matrix \f$\mathbf{Q} = \left[\tilde\Sigma_{jk}\right]\f$ with entries,
!>  \f{equation}{
!>      \tilde\Sigma_{jk} = \frac{1}{n} \sum_{i=1}^{n} \left( x_{ij} - \hat\mu_j \right) \left( x_{ik} - \hat\mu_k \right) ~,
!>  \f}
!>  where \f$n\f$ is the number of observations in the sample, \f$\hat\mu\f$ is the sample mean vector, and \f$\Sigma_{jk}\f$ is an
!>  estimate of the covariance between the \f$j\f$th variable and the kth variable of the population underlying the data.<br>
!>
!>  The diagonal elements of the matrix \f$\tilde\Sigma_{jj}\f$ are known as the **sample variance**.
!>
!>  Biased sample covariance
!>  ------------------------
!>
!>  The above formula yields a biased estimate of the covariance matrix of the sample.<br>
!>  Intuitively, the sample covariance relies on the difference between each observation and the sample mean,
!>  but the sample mean is slightly correlated with each observation since it is defined in terms of all observations.<br>
!>  Therefore, unless the sample mean is known a priori, the above equation yields a biased estimate of the covariance
!>  with sample mean as a proxy for the true mean of the population.<br>
!>  Note that the bias is noticeable only when the sample size is small (e.g., \f$<10\f$).<br>
!>
!>  Unbiased sample covariance
!>  --------------------------
!>
!>  A popular fix to the definition of sample covariance to remove its bias is to apply the
!>  [Bessel correction](https://en.wikipedia.org/wiki/Bessel%27s_correction) to the equation above,
!>  yielding the unbiased covariance estimate as,
!>  \f{eqnarray}{
!>      \hat\Sigma_{jk}
!>      &=& \frac{\xi}{n} \sum_{i=1}^{n} \left( x_{ij} - \hat\mu_j \right) \left( x_{ik} - \hat\mu_k \right) ~,
!>      &=& \frac{1}{n - 1} \sum_{i=1}^{n} \left( x_{ij} - \hat\mu_j \right) \left( x_{ik} - \hat\mu_k \right) ~,
!>  \f}
!>  where \f$\xi = \frac{n}{n - 1}\f$ is the **Bessel bias correction factor**.<br>
!>
!>  Biased weighted sample covariance
!>  ---------------------------------
!>
!>  \f{equation}{
!>      \tilde{\Sigma}^w_{jk} = \frac{ \sum_{i = 1}^{n} \left( x_{ij} - \hat\mu^w_j \right) \left( x_{ik} - \hat\mu^w_k \right) } { \left( \sum_{i=1}^{n} w_i \right) } ~.
!>  \f}
!>  where `n = nsam` is the number of observations in the sample, \f$w_i\f$ are the weights of individual data points,
!>  the superscript \f$^w\f$ signifies the sample weights, and \f$\hat\mu^w\f$ is the **weighted mean** of the sample.<br>
!>  When the sample size is small, the above equation yields a biased estimate of the covariance.<br>
!>
!>  **Unbiased weighted sample covariance**<br>
!>
!>  There is no unique generic equation for the **unbiased covariance** of a **weighted sample**.<br>
!>  However, depending on the types of the weights involved, a few popular definitions exist.<br>
!>
!>  <ol>
!>      <li>    The **unbiased covariance** of a sample with **frequency**, **count**, or **repeat** weights can be computed via the following equation,
!>              \f{equation}{
!>                  \hat\Sigma^w_{jk} = \frac{ \sum_{i = 1}^{n} \left( x_{ij} - \hat\mu^w_j \right) \left( x_{ik} - \hat\mu^w_k \right) } { \left( \sum_{i=1}^{n} w_i \right) - 1} ~.
!>              \f}
!>              [Frequency weights](@ref pm_sampleWeight::fweight_type) represent the number of duplications of each observation in the sample whose population covariance is to be estimated.<br>
!>              Therefore, the **frequency weights** are expected to be **integers** or **whole numbers**.<br>
!>      <li>    The **unbiased covariance** of a sample with **reliability weights**, also sometimes confusingly known as **probability weights** or **importance weights**,
!>              can be computed by the following equation,
!>              \f{equation}{
!>                  \hat\Sigma^w_{jk} = \frac{ \sum_{i=1}^{n} w_i } { \left( \sum_{i=1}^{n} w_i \right)^2 - \left( \sum_{i=1}^{n} w_i^2 \right) } \sum_{i = 1}^{n} \left( x_{ij} - \hat\mu^w_j \right) \left( x_{ik} - \hat\mu^w_k \right) ~.
!>              \f}
!>              <ol>
!>                  <li>    [Reliability weights](@ref pm_sampleWeight::rweight_type) weights, also known as **reliability weights** or **sampling weights**
!>                          represent the probability of a case (or subject) being selected into the sample from a population.<br>
!>                  <li>    Application of the term *unbiased* to the above equation is controversial as some believe that bias
!>                          cannot be correct without the knowledge of the sample size, which is lost in normalized weights.<br>
!>                  <li>    Reliability weights are frequently (but not necessarily) normalized, meaning that \f$\sum^{i = 1}_{n} w_i = 1\f$.<br>
!>              </ol>
!>  </ol>
!>
!>  Covariance matrix vs. correlation matrix
!>  ----------------------------------------
!>
!>  The covariance matrix \f$\Sigma\f$ is related to the [correlation matrix](@ref pm_sampleCor) \f$\rho\f$ by the following equation,
!>  \f{equation}{
!>      \Sigma_{ij} = \rho_{ij} \times \sigma_{i} \times \sigma_{j} ~,
!>  \f}
!>  where \f$\Sigma\f$ represents the covariance matrix, \f$\rho\f$ represents the correlation matrix, and \f$\sigma\f$ represents the standard deviations.<br>
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
!>  Box and Tiao, 1973, *Bayesian Inference in Statistical Analysis*, Page 421.<br>
!>  *Updating mean and variance estimates: an improved method*, D.H.D. West, 1979.<br>
!>  Geisser and Cornfield, 1963, *Posterior distributions for multivariate normal parameters*.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{setCov_vs_setCovMean, The runtime performance of [setCov](@ref pm_sampleCov::setCov) vs. [setCovMean](@ref pm_sampleCov::setCovMean).}
!>  \include{lineno} benchmark/pm_sampleCov/setCov_vs_setCovMean/main.F90
!>  \compilefb{setCov_vs_setCovMean}
!>  \postprocb{setCov_vs_setCovMean}
!>  \include{lineno} benchmark/pm_sampleCov/setCov_vs_setCovMean/main.py
!>  \visb{setCov_vs_setCovMean}
!>  \image html benchmark/pm_sampleCov/setCov_vs_setCovMean/benchmark.setCov_vs_setCovMean.runtime.png width=1000
!>  \image html benchmark/pm_sampleCov/setCov_vs_setCovMean/benchmark.setCov_vs_setCovMean.runtime.ratio.png width=1000
!>  \moralb{setCov_vs_setCovMean}
!>      -#  The procedures under the generic interface [setCov](@ref pm_sampleCov::setCov) take the sample mean as input and return the covariance matrix.<br>
!>      -#  The procedures under the generic interface [setCovMean](@ref pm_sampleCov::setCovMean) compute both the sample mean and covariance matrix in one pass.<br>
!>      -#  The performance of the two methods appears to depend significantly on the compiler used.<br>
!>      -#  But in general, the one-pass algorithm of [setCovMean](@ref pm_sampleCov::setCovMean) appears to
!>          perform equally or slightly better than the two-pass algorithm of [setCov](@ref pm_sampleCov::setCov).<br>
!>
!>  \test
!>  [test_pm_sampleCov](@ref test_pm_sampleCov)
!>
!>  \bug
!>  \status See \unresolved, See [this page](https://groups.google.com/g/comp.lang.fortran/c/NDE6JKTFbNU) for more information.<br>
!>  \source \gfortran
!>  \desc
!>  Ideally, there should be only one generic interface in this module for computing the biased/corrected/weighted variance.<br>
!>  This requires ability to resolve the different weight types, which requires custom derived types for weights.<br>
!>  Fortran PDTs are ideal for such use cases. However, the implementation of PDTs is far from complete in \gfortran.<br>
!>  \remedy
!>  Given that the importance of \gfortran support, separate generic interfaces were instead developed for different sample weight types.<br>
!>  Once the \gfortran PDT bugs are resolved, the [getVar](@ref pm_sampleVar::getVar) generic interface can be extended to serve as a
!>  high-level wrapper for the weight-specific generic interfaces in this module.<br>
!>
!>  \todo
!>  \pmed
!>  The inclusion of bias correction in the calculation of covariance is a
!>  frequentist abomination and shenanigan that must be eliminated in the future.<br>
!>  The correction factor should be computed separately from the actual covariance calculation.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Nov 24, 2020, 4:19 AM, Dallas, TX<br>
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX<br>
!>  \AmirShahmoradi, Monday March 6, 2017, 2:48 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleCov

    use pm_kind, only: SK, IK, LK, RK
    use pm_sampleWeight, only: weight_type
    use pm_sampleWeight, only: fweight_type, fweight
    use pm_sampleWeight, only: rweight_type, rweight
    use pm_matrixSubset, only: uppDia_type, uppDia
    use pm_matrixSubset, only: lowDia_type, lowDia
    use pm_sampleVar, only: getVarCorrection

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleCov"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (optionally unbiased) covariance matrix of a pair of (potentially weighted) time series `x(1:nsam)` and `y(1:nsam)` or
    !>  of an input (potentially weighted) array of shape `(ndim, nsam)` or `(nsam, ndim)` where `ndim` is the number of data dimensions
    !>  (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  This generic interface performs one of the following computational tasks:<br>
    !>  <ol>
    !>      <li>    compute the covariance matrix corresponding to an input correlation matrix and vector of standard deviations in arbitrary `ndim` dimensions.<br>
    !>      <li>    compute the sample covariance matrix of a random sample of `nsam` observations in arbitrary `ndim` dimensional space.<br>
    !>      <li>    compute the sample covariance matrix of a pair of `nsam` observations in two separated data vectors `x` and `y`.<br>
    !>  </ol>
    !>  See the documentation of the parent module [pm_sampleCov](@ref pm_sampleCov) for algorithmic details and sample covariance matrix definition.<br>
    !>
    !   \param[in]  subset      :   The input scalar constant argument that can be any of the following:<br>
    !                               <ol>
    !                                   <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the output covariance matrix must be computed.<br>
    !                                   <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the output covariance matrix must be computed.<br>
    !                               </ol>
    !                               This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !                               (**optional**. If missing, the full covariance matrix will be output. It can be present **if and only if** the input arguments `x` and `y` are missing.)
    !>  \param[in]  cor         :   The input positive semi-definite square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the input `cov`,
    !>                              representing the correlation matrix based upon which the output covariance matrix `cov` is to be computed.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input arguments `std` and `subsetr` are present and the rest are missing.)
    !>  \param[in]  subsetr     :   The input scalar constant argument that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input correlation matrix must be used.<br>
    !>                                  <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input correlation matrix must be used.<br>
    !>                              </ol>
    !>                              This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                              Although the allowed subset constants imply the use of the diagonal elements, they are by definition assumed to be `1` and therefore are not referenced in the algorithm.
    !>                              (**optional**. If missing, only the upper triangle will be used. It **must** be present **if and only if** the input arguments `cor` is present.)
    !>  \param[in]  std         :   The input positive vector of shape `(1:ndim)` of type `real` of the same kind as the input `cov`,
    !>                              containing the standard deviations of the `ndim` data attributes based upon which the output covariance matrix `cov` is to be computed.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `cor` is present.)
    !>  \param[in]  x           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                              containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `y` is present and `sample` is missing.)
    !>  \param[in]  y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                              containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `x` is present and `sample` is missing.)
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `cov`,
    !>                              containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                              If `sample` is a matrix, then the direction along which the covariance matrix is computed is dictated by the input argument `dim`.<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `dim` is present and `x` and `y` are missing.)
    !>  \param[in]  dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the covariance matrix must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input argument `sample` is present.)
    !>  \param[in]  weight      :   The `contiguous` vector of length `nsam` of,
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK, or
    !>                                  <li>    type `real` of the same kind as the output `cov`,
    !>                              </ol>
    !>                              containing the corresponding weights of individual `nsam` observations in `sample`.<br>
    !>                              (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled). It can be present **only if** the input arguments `cor`, `subsetr`, and `std` are missing.)
    !>  \param[in]  correction  :   The input scalar object that can be the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [fweight](@ref pm_sampleWeight::fweight) or an object of type [fweight](@ref pm_sampleWeight::fweight)
    !>                                          implying a bias correction based on the assumption of **frequency weights** for the sample observations, even if the `weight` argument is missing.<br>
    !>                                          This is the most popular type of correction, also known as the [Bessel correction](https://en.wikipedia.org/wiki/Bessel%27s_correction).<br>
    !>                                  <li>    The constant [rweight](@ref pm_sampleWeight::rweight) or an object of type [rweight_type](@ref pm_sampleWeight::rweight_type)
    !>                                          implying a bias correction based on the assumption of **reliability weights** for the sample observations.<br>
    !>                              </ol>
    !>                              (**optional**. If missing, no bias-correction will be applied to the output `cov`.<br>
    !>                              It can be present **only if** the input arguments `cor`, `subsetr`, and `std` are missing.)
    !>
    !>  \return
    !>  `cov`                   :   The output **positive semi-definite** square matrix of shape `(1:ndim, 1:ndim)` of,
    !>                              <ol>
    !>                                  <li>     type `complex` of kind \CKALL,
    !>                                  <li>     type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the computed covariance matrix.<br>
    !>                              When the input arguments `x` and `y` are present, then `ndim == 2` holds by definition and the output covariance is of the form,
    !>                              \f{equation}{
    !>                                  \ms{cov} =
    !>                                  \begin{bmatrix}
    !>                                      \sigma_{xx} && \sigma_{yx} \\
    !>                                      \sigma_{xy} && \sigma_{yy}
    !>                                  \end{bmatrix}
    !>                              \f}
    !>
    !>  \interface{getCov}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: getCov
    !>
    !>      ! correlation matrix to covariance matrix.
    !>
    !>      cov(1:ndim, 1:ndim) = getCov(cor(1:ndim, 1:ndim), subsetr, std(1:ndim))
    !>
    !>      ! XY time series covariance matrix.
    !>
    !>      cov(1:2, 1:2) = getCov(x(1:nsam), y(1:nsam)                , correction = correction) ! full covariance matrix.
    !>      cov(1:2, 1:2) = getCov(x(1:nsam), y(1:nsam), weight(1:nsam), correction = correction) ! full covariance matrix.
    !>
    !>      ! sample covariance matrix.
    !>
    !>      cov(1:ndim, 1:ndim) = getCov(sample(:,:), dim                , correction = correction) ! full covariance matrix.
    !>      cov(1:ndim, 1:ndim) = getCov(sample(:,:), dim, weight(1:nsam), correction = correction) ! full covariance matrix.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  This generic interface is merely a flexible wrapper around the generic `subroutine` interface [setCov](@ref pm_sampleCov::setCov).<br>
    !>  As such, all conditions and warnings associated with [setCov](@ref pm_sampleCov::setCov) equally hold for this generic interface.<br>
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  Note the effects of bias-correction in computing the covariance become noticeable only for very sample sample sizes (i.e., when `nsam` is small).<br>
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{getCov}
    !>  \include{lineno} example/pm_sampleCov/getCov/main.F90
    !>  \compilef{getCov}
    !>  \output{getCov}
    !>  \include{lineno} example/pm_sampleCov/getCov/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.8}
    !>  \desc
    !>  Intel ifort (possibly only with heap memory allocation) yields a segfault error in OpenMP parallel code at `pm_sampleCov@routines.inc.F90:238`.<br>
    !>  The root cause of the segmentation fault was determined to be due to the use of `do concurrent` construct in the implementation.<br>
    !>  \remedy
    !>  For now, all do concurrent statements are converted back to normal do loops.<br>
    !>
    !>  \final{getCov}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! cor2cov.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovCorStd_ULD_UXD_CK5(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovCorStd_ULD_UXD_CK4(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovCorStd_ULD_UXD_CK3(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovCorStd_ULD_UXD_CK2(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovCorStd_ULD_UXD_CK1(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovCorStd_ULD_UXD_RK5(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovCorStd_ULD_UXD_RK4(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovCorStd_ULD_UXD_RK3(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovCorStd_ULD_UXD_RK2(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovCorStd_ULD_UXD_RK1(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovCorStd_ULD_XLD_CK5(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovCorStd_ULD_XLD_CK4(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovCorStd_ULD_XLD_CK3(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovCorStd_ULD_XLD_CK2(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovCorStd_ULD_XLD_CK1(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)                                                :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovCorStd_ULD_XLD_RK5(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovCorStd_ULD_XLD_RK4(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovCorStd_ULD_XLD_RK3(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovCorStd_ULD_XLD_RK2(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovCorStd_ULD_XLD_RK1(cor, subsetr, std) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovCorStd_ULD_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)                                                   :: cov(size(cor, 1, IK), size(cor,2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! no weight: XY.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovWNO_XY_CK5(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovWNO_XY_CK4(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovWNO_XY_CK3(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovWNO_XY_CK2(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovWNO_XY_CK1(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovWNO_XY_RK5(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovWNO_XY_RK4(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovWNO_XY_RK3(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovWNO_XY_RK2(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovWNO_XY_RK1(x, y, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCov

    ! integer weight: XY.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovWTI_XY_CK5(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovWTI_XY_CK4(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovWTI_XY_CK3(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovWTI_XY_CK2(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovWTI_XY_CK1(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovWTI_XY_RK5(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovWTI_XY_RK4(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovWTI_XY_RK3(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovWTI_XY_RK2(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovWTI_XY_RK1(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCov

    ! real weight: XY.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovWTR_XY_CK5(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovWTR_XY_CK4(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovWTR_XY_CK3(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovWTR_XY_CK2(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovWTR_XY_CK1(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: x(:), y(:)
        complex(TKG)                                                    :: cov(2, 2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovWTR_XY_RK5(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovWTR_XY_RK4(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovWTR_XY_RK3(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovWTR_XY_RK2(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovWTR_XY_RK1(x, y, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: x(:), y(:)
        real(TKG)                                                       :: cov(2, 2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! no weight: ULD.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovWNO_ULD_CK5(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovWNO_ULD_CK4(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovWNO_ULD_CK3(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovWNO_ULD_CK2(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovWNO_ULD_CK1(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovWNO_ULD_RK5(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovWNO_ULD_RK4(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovWNO_ULD_RK3(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovWNO_ULD_RK2(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovWNO_ULD_RK1(sample, dim, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWNO_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCov

    ! integer weight: ULD.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovWTI_ULD_CK5(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovWTI_ULD_CK4(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovWTI_ULD_CK3(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovWTI_ULD_CK2(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovWTI_ULD_CK1(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovWTI_ULD_RK5(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovWTI_ULD_RK4(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovWTI_ULD_RK3(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovWTI_ULD_RK2(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovWTI_ULD_RK1(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTI_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        integer(IK)             , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCov

    ! real weight: ULD.

    interface getCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovWTR_ULD_CK5(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovWTR_ULD_CK4(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovWTR_ULD_CK3(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovWTR_ULD_CK2(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovWTR_ULD_CK1(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        complex(TKG)            , intent(in)    , contiguous            :: sample(:,:)
        complex(TKG)                                                    :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovWTR_ULD_RK5(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovWTR_ULD_RK4(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovWTR_ULD_RK3(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovWTR_ULD_RK2(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovWTR_ULD_RK1(sample, dim, weight, correction) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovWTR_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                            :: dim
        class(weight_type)      , intent(in)                , optional  :: correction
        real(TKG)               , intent(in)    , contiguous            :: weight(:)
        real(TKG)               , intent(in)    , contiguous            :: sample(:,:)
        real(TKG)                                                       :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCov

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the covariance matrix corresponding to the input (potentially weighted) correlation matrix or
    !>  return the **biased** sample covariance matrix of the input array of shape `(ndim, nsam)` or `(nsam, ndim)` or
    !>  a pair of (potentially weighted) time series `x(1:nsam)` and `y(1:nsam)` where `ndim` is the number of data dimensions
    !>  (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  This generic interface performs either of the following two computational tasks:<br>
    !>  <ol>
    !>      <li>    compute the covariance matrix corresponding to an input correlation matrix and vector of standard deviations in arbitrary `ndim` dimensions.<br>
    !>      <li>    compute the sample covariance matrix of a random sample of `nsam` observations in arbitrary `ndim` dimensional space.<br>
    !>  </ol>
    !>  See the documentation of the parent module [pm_sampleCov](@ref pm_sampleCov) for algorithmic details and sample covariance matrix definition.<br>
    !>
    !>  \param[inout]   cov         :   The output positive semi-definite square matrix of shape `(1:ndim, 1:ndim)` of,
    !>                                  <ol>
    !>                                      <li>     type `complex` of kind \CKALL,
    !>                                      <li>     type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  representing the covariance matrix corresponding to the input `sample` or `cor`, whichever is present.<br>
    !>                                  On output,<br>
    !>                                  <ol>
    !>                                      <li>    If the input arguments `x` and `y` are present, then `ndim == 2` holds by definition and the output covariance is of the form
    !>                                              \f{equation}{
    !>                                                  \ms{cov} =
    !>                                                  \begin{bmatrix}
    !>                                                      \sigma_{xx} && \sigma_{yx} \\
    !>                                                      \sigma_{xy} && \sigma_{yy}
    !>                                                  \end{bmatrix}
    !>                                              \f}
    !>                                      <li>    Otherwise, is `subset` is present, then only the specified input `subset` will be overwritten with the covariance matrix.<br>
    !>                                              Any elements not in the specified input `subset` remain intact.<br>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar constant argument that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the output covariance matrix must be computed.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the output covariance matrix must be computed.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `x` and `y` are both missing.)
    !>  \param[in]      cor         :   The input positive semi-definite square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the input `cov`,
    !>                                  representing the correlation matrix based upon which the output covariance matrix `cov` is to be computed.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `std` and `subsetr` are present and the rest are missing.)
    !>  \param[in]      subsetr     :   The input scalar constant argument that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input correlation matrix must be used.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input correlation matrix must be used.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                                  Although the allowed subset constants imply the use of the diagonal elements, they are by definition assumed to be `1` and therefore are not referenced in the algorithm.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `subset`, `cor`, and `std` are present and the rest are missing.)
    !>  \param[in]      std         :   The input positive vector of shape `(1:ndim)` of type `real` of the same kind as the input `cov`,
    !>                                  containing the standard deviations of the `ndim` data attributes based upon which the output covariance matrix `cov` is to be computed.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `cor` and `subsetr` are present and the rest are missing.)
    !>  \param[in]      mean        :   The input `contiguous` vector of shape `(ndim)` of the same type and kind as the output `cov` containing the `sample` mean.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `sample` is missing and `x` and `y` are present, then `mean` must be a vector of size `ndim = 2`.
    !>                                      <li>    If the input `sample` is present and is a 2D array, then `mean` must be a vector of size `ndim = size(sample, 3 - dim)`
    !>                                              (i.e., computed along the specified input dimension `dim`).<br>
    !>                                  </ol>
    !>                                  The mean vector can be readily computed via [getMean](@ref pm_sampleMean::getMean) or [setMean](@ref pm_sampleMean::setMean).<br>
    !>                                  (**optional**. default = [getFilled(0., ndim)](@ref pm_arrayFill::getFilled) or `[0., 0.]` depending on the rank of the input `sample` and `dim` or the presence of `x` and `y`.)
    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `cor` and `sample` are missing.)
    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `cor` and `sample` are missing.)
    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `cov`,
    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                                  If `sample` is a matrix, then the direction along which the covariance matrix is computed is dictated by the input argument `dim`.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `dim` is present and `x` and `y` are missing.)
    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the covariance matrix must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present.)
    !>  \param[in]      weight      :   The `contiguous` vector of length `nsam` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, or
    !>                                      <li>    type `real` of the same kind as the kind of the output `cov`,
    !>                                  </ol>
    !>                                  containing the corresponding weights of individual `nsam` observations in `sample`.<br>
    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled). It can be present **if and only if** the input arguments `sample` or `x` and `y` are present.)
    !>  \param[in]      weisum      :   The input scalar of the same type and kind as the input `weight` containing `sum(weight)`.<br>
    !>                                  This quantity is a byproduct of computing the mean of a sample and is automatically returned by [setMean](@ref pm_sampleMean::setMean).<br>
    !>                                  (**optional**. It **must** be present **if and only if** the `weight` argument is present **and** the input arguments `cor`, `subsetr`, and `std` are missing.)
    !>
    !>  \interface{setCov}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: setCov, lowDia, uppDia
    !>
    !>      ! correlation matrix to covariance matrix.
    !>
    !>      call setCov(cov(1:ndim, 1:ndim), subset, cor(1:ndim, 1:ndim), subsetr, std(1:ndim))
    !>
    !>      ! XY time series covariance matrix.
    !>
    !>      call setCov(cov(1:2, 1:2)           , x(1:nsam), y(1:nsam))
    !>      call setCov(cov(1:2, 1:2), mean(1:2), x(1:nsam), y(1:nsam))
    !>      call setCov(cov(1:2, 1:2)           , x(1:nsam), y(1:nsam), weight(1:nsam))
    !>      call setCov(cov(1:2, 1:2), mean(1:2), x(1:nsam), y(1:nsam), weight(1:nsam), weisum)
    !>
    !>      ! sample covariance matrix.
    !>
    !>      call setCov(cov(1:ndim, 1:ndim), subset              , sample(:,:), dim)
    !>      call setCov(cov(1:ndim, 1:ndim), subset, mean(1:ndim), sample(:,:), dim)
    !>      call setCov(cov(1:ndim, 1:ndim), subset,             , sample(:,:), dim, weight(1:nsam))
    !>      call setCov(cov(1:ndim, 1:ndim), subset, mean(1:ndim), sample(:,:), dim, weight(1:nsam), weisum)
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
    !>  \note
    !>  Note the effects of bias-correction in computing the covariance matrix become noticeable only for very sample sample sizes (i.e., when `nsam` is small).<br>
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{setCov}
    !>  \include{lineno} example/pm_sampleCov/setCov/main.F90
    !>  \compilef{setCov}
    !>  \output{setCov}
    !>  \include{lineno} example/pm_sampleCov/setCov/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setCov_dim1_vs_dim2, The runtime performance of [setCov](@ref pm_sampleCov::setCov) along different sample dimensions.}
    !>  \include{lineno} benchmark/pm_sampleCov/setCov_dim1_vs_dim2/main.F90
    !>  \compilefb{setCov_dim1_vs_dim2}
    !>  \postprocb{setCov_dim1_vs_dim2}
    !>  \include{lineno} benchmark/pm_sampleCov/setCov_dim1_vs_dim2/main.py
    !>  \visb{setCov_dim1_vs_dim2}
    !>  \image html benchmark/pm_sampleCov/setCov_dim1_vs_dim2/benchmark.setCov_dim1_vs_dim2.runtime.png width=1000
    !>  \image html benchmark/pm_sampleCov/setCov_dim1_vs_dim2/benchmark.setCov_dim1_vs_dim2.runtime.ratio.png width=1000
    !>  \moralb{setCov_dim1_vs_dim2}
    !>      -#  The procedures under the generic interface [setCov](@ref pm_sampleCov::setCov) can compute the covariance under two different sample axes.<br>
    !>      -#  Recall that C is a row-major language while Fortran is a column-major language.<br>
    !>      -#  As such, one would expect the computations for a sample whose observations are stored along the second axis would be faster in the Fortran programming language.<br>
    !>      -#  However, such an expectation does not appear to hold at all times and appears to depend significantly on the computing architecture and the number of data attributes involved.<br>
    !>      -#  The higher the number of data attributes, the more likely the computations along the second axis of `sample` will be faster.<br>
    !>      -#  Note that for small number of data attributes, the computations along the second data axis involve a small
    !>          loop that has significant computational cost due to the implicit branching involved in the loop.<br>
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \todo
    !>  \pmed
    !>  The examples of this generic interface should be extended to corrected weighted covariance matrices.<br>
    !>
    !>  \final{setCov}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX

    ! cor2cov.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_CK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_CK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_CK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_CK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_CK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_RK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_RK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_RK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_RK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovCorStd_UXD_UXD_RK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_CK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_CK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_CK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_CK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_CK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_RK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_RK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_RK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_RK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovCorStd_UXD_XLD_RK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_UXD_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(uppDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_CK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_CK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_CK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_CK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_CK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_RK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_RK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_RK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_RK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovCorStd_XLD_UXD_RK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(uppDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_CK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_CK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_CK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_CK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_CK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        complex(TKG)        , intent(in)    , contiguous            :: cor(:,:)
        complex(TKG)        , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_RK5(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_RK4(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_RK3(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_RK2(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovCorStd_XLD_XLD_RK1(cov, subset, cor, subsetr, std)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovCorStd_XLD_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(lowDia_type)   , intent(in)                            :: subset
        type(lowDia_type)   , intent(in)                            :: subsetr
        real(TKG)           , intent(in)    , contiguous            :: std(:)
        real(TKG)           , intent(in)    , contiguous            :: cor(:,:)
        real(TKG)           , intent(inout) , contiguous            :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! XY - no weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWNO_XY_CK5(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWNO_XY_CK4(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWNO_XY_CK3(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWNO_XY_CK2(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWNO_XY_CK1(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWNO_XY_RK5(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWNO_XY_RK4(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWNO_XY_RK3(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWNO_XY_RK2(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWNO_XY_RK1(cov, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWNO_XY_CK5(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWNO_XY_CK4(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWNO_XY_CK3(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWNO_XY_CK2(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWNO_XY_CK1(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWNO_XY_RK5(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWNO_XY_RK4(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWNO_XY_RK3(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWNO_XY_RK2(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWNO_XY_RK1(cov, mean, x, y)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    ! XY - integer weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWTI_XY_CK5(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWTI_XY_CK4(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWTI_XY_CK3(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWTI_XY_CK2(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWTI_XY_CK1(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWTI_XY_RK5(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWTI_XY_RK4(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWTI_XY_RK3(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWTI_XY_RK2(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWTI_XY_RK1(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWTI_XY_CK5(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWTI_XY_CK4(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWTI_XY_CK3(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWTI_XY_CK2(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWTI_XY_CK1(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWTI_XY_RK5(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWTI_XY_RK4(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWTI_XY_RK3(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWTI_XY_RK2(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWTI_XY_RK1(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    ! XY - real weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWTR_XY_CK5(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWTR_XY_CK4(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWTR_XY_CK3(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWTR_XY_CK2(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWTR_XY_CK1(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWTR_XY_RK5(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWTR_XY_RK4(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWTR_XY_RK3(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWTR_XY_RK2(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWTR_XY_RK1(cov, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWTR_XY_CK5(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWTR_XY_CK4(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWTR_XY_CK3(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWTR_XY_CK2(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWTR_XY_CK1(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWTR_XY_RK5(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWTR_XY_RK4(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWTR_XY_RK3(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWTR_XY_RK2(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWTR_XY_RK1(cov, mean, x, y, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! UXD - no weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_CK5(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_CK4(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_CK3(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_CK2(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_CK1(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_RK5(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_RK4(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_RK3(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_RK2(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWNO_UXD_RK1(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_CK5(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_CK4(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_CK3(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_CK2(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_CK1(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_RK5(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_RK4(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_RK3(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_RK2(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWNO_UXD_RK1(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    ! UXD - integer weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_CK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_CK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_CK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_CK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_CK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_RK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_RK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_RK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_RK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWTI_UXD_RK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_CK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_CK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_CK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_CK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_CK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_RK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_RK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_RK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_RK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWTI_UXD_RK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    ! UXD - real weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_CK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_CK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_CK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_CK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_CK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_RK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_RK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_RK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_RK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWTR_UXD_RK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_CK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_CK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_CK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_CK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_CK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_RK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_RK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_RK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_RK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWTR_UXD_RK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! XLD - no weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_CK5(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_CK4(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_CK3(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_CK2(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_CK1(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_RK5(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_RK4(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_RK3(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_RK2(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWNO_XLD_RK1(cov, subset, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWNO_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_CK5(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_CK4(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_CK3(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_CK2(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_CK1(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_RK5(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_RK4(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_RK3(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_RK2(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWNO_XLD_RK1(cov, subset, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWNO_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    ! XLD - integer weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_CK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_CK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_CK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_CK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_CK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_RK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_RK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_RK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_RK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWTI_XLD_RK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTI_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_CK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_CK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_CK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_CK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_CK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_RK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_RK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_RK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_RK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWTI_XLD_RK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTI_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

    ! XLD - real weight.

    interface setCov

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_CK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_CK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_CK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_CK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_CK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_RK5(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_RK4(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_RK3(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_RK2(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovOrgWTR_XLD_RK1(cov, subset, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovOrgWTR_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_CK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_CK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_CK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_CK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_CK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        complex(TKG)            , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_RK5(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_RK4(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_RK3(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_RK2(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovAvgWTR_XLD_RK1(cov, subset, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovAvgWTR_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)                        :: weisum
        real(TKG)               , intent(in)    , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCov

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the covariance matrix and mean vector corresponding to the input (potentially weighted) input `sample` of shape `(ndim, nsam)` or `(nsam, ndim)`
    !>  or a pair of (potentially weighted) time series `x(1:nsam)` and `y(1:nsam)` where `ndim` is the number of data dimensions (the number of data attributes)
    !>  and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_sampleCov](@ref pm_sampleCov) for algorithmic details and sample covariance matrix definition.<br>
    !>
    !>  \param[inout]   cov         :   The output positive semi-definite square matrix of shape `(1:ndim, 1:ndim)` of,
    !>                                  <ol>
    !>                                      <li>     type `complex` of kind \CKALL,
    !>                                      <li>     type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  representing the covariance matrix corresponding to the input `sample` or the pair `x` and `y`, whichever is present.<br>
    !>                                  On output,<br>
    !>                                  <ol>
    !>                                      <li>    If the input arguments `x` and `y` are present, then `ndim == 2` holds by definition and the output covariance is of the form
    !>                                              \f{equation}{
    !>                                                  \ms{cov} =
    !>                                                  \begin{bmatrix}
    !>                                                      \sigma_{xx} && \sigma_{yx} \\
    !>                                                      \sigma_{xy} && \sigma_{yy}
    !>                                                  \end{bmatrix}
    !>                                              \f}
    !>                                      <li>    Otherwise, is `subset` is present, then only the specified input `subset` will be overwritten with the covariance matrix.<br>
    !>                                              Any elements not in the specified input `subset` remain intact.<br>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar constant argument that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the output covariance matrix must be computed.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the output covariance matrix must be computed.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input arguments `x` and `y` are both missing.)
    !>  \param[out]     mean        :   The output `contiguous` vector of shape `(ndim)` of the same type and kind as the output `cov` containing the `sample` mean.<br>
    !>                                  <ol>
    !>                                      <li>    If the input `sample` is missing and `x` and `y` are present, then `mean` must be a vector of size `ndim = 2`.
    !>                                      <li>    If the input `sample` is present and is a 2D array, then `mean` must be a vector of size `ndim = size(sample, 3 - dim)`
    !>                                              (i.e., computed along the specified input dimension `dim`).<br>
    !>                                  </ol>
    !>  \param[in]      x           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                                  containing the first attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `y` is present and `sample` is missing.)
    !>  \param[in]      y           :   The input `contiguous` vector of shape `(nsam)` of the same type and kind as the output `cov`,
    !>                                  containing the second attribute `x` of the observational sample, where `nsam` is the number of observations in the sample.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `x` is present and `sample` is missing.)
    !>  \param[in]      sample      :   The input `contiguous` array of shape `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `cov`,
    !>                                  containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                                  If `sample` is a matrix, then the direction along which the covariance matrix is computed is dictated by the input argument `dim`.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the input argument `dim` is present and `x` and `y` are missing.)
    !>  \param[in]      dim         :   The input scalar `integer` of default kind \IK indicating the dimension of `sample` along which the covariance matrix must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `dim = 1`, the input `sample` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                      <li>    If `dim = 2`, the input `sample` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                                  </ol>
    !>                                  (**optional**. It must be present **if and only if** the input argument `sample` is present.)
    !>  \param[in]      weight      :   The `contiguous` vector of length `nsam` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK, or
    !>                                      <li>    type `real` of the same kind as the kind of the output `cov`,
    !>                                  </ol>
    !>                                  containing the corresponding weights of individual `nsam` observations in `sample`.<br>
    !>                                  (**optional**. default = [getFilled(1, nsam)](@ref pm_arrayFill::getFilled). It can be present **if and only if** the input arguments `sample` or `x` and `y` are present.)
    !>  \param[out]     weisum      :   The output scalar of the same type and kind as the input `weight`, containing `sum(weight)`.<br>
    !>                                  (**optional**. It **must** be present **if and only if** the `weight` argument is present.)
    !>  \param[in]      meang        :  The input vector of the same type, kind, and rank as the output `cov` of shape `(1:ndim)` containing the best guess for the `sample` mean.<br>
    !>                                  If no good guess is known a priori, `meang` can be set to any observation in `sample`.<br>
    !>                                  See the example usage below.<br>
    !>
    !>  \interface{setCovMean}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: setCovMean
    !>
    !>      ! XY time series covariance matrix.
    !>
    !>      call setCovMean(cov(1:2, 1:2), mean(1:2), x(1:nsam), y(1:nsam)                        , meang(1:ndim))
    !>      call setCovMean(cov(1:2, 1:2), mean(1:2), x(1:nsam), y(1:nsam), weight(1:nsam), weisum, meang(1:ndim))
    !>
    !>      ! sample covariance matrix.
    !>
    !>      call setCovMean(cov(1:ndim, 1:ndim), subset, mean(1:ndim)   , sample(:,:), dim                          , meang(1:ndim))
    !>      call setCovMean(cov(1:ndim, 1:ndim), subset, mean(1:ndim)   , sample(:,:), dim, weight(1:nsam), weisum  , meang(1:ndim))
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
    !>  The condition `1 < size(x)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 < size(y)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  Note the effects of bias-correction in computing the covariance matrix become noticeable only for very sample sample sizes (i.e., when `nsam` is small).<br>
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{setCovMean}
    !>  \include{lineno} example/pm_sampleCov/setCovMean/main.F90
    !>  \compilef{setCovMean}
    !>  \output{setCovMean}
    !>  \include{lineno} example/pm_sampleCov/setCovMean/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setCovMean_dim1_vs_dim2, The runtime performance of [setCovMean](@ref pm_sampleCov::setCovMean) along different sample dimensions.}
    !>  \include{lineno} benchmark/pm_sampleCov/setCovMean_dim1_vs_dim2/main.F90
    !>  \compilefb{setCovMean_dim1_vs_dim2}
    !>  \postprocb{setCovMean_dim1_vs_dim2}
    !>  \include{lineno} benchmark/pm_sampleCov/setCovMean_dim1_vs_dim2/main.py
    !>  \visb{setCovMean_dim1_vs_dim2}
    !>  \image html benchmark/pm_sampleCov/setCovMean_dim1_vs_dim2/benchmark.setCovMean_dim1_vs_dim2.runtime.png width=1000
    !>  \image html benchmark/pm_sampleCov/setCovMean_dim1_vs_dim2/benchmark.setCovMean_dim1_vs_dim2.runtime.ratio.png width=1000
    !>  \moralb{setCovMean_dim1_vs_dim2}
    !>      -#  The procedures under the generic interface [setCovMean](@ref pm_sampleCov::setCovMean) can compute the covariance under two different sample axes.<br>
    !>      -#  Recall that C is a row-major language while Fortran is a column-major language.<br>
    !>      -#  As such, one would expect the computations for a sample whose observations are stored along the second axis would be faster in the Fortran programming language.<br>
    !>      -#  However, such an expectation does not appear to hold at all times and appears to depend significantly on the computing architecture and the number of data attributes involved.<br>
    !>      -#  The higher the number of data attributes, the more likely the computations along the second axis of `sample` will be faster.<br>
    !>      -#  Note that for small number of data attributes, the computations along the second data axis involve a small
    !>          loop that has significant computational cost due to the implicit branching involved in the loop.<br>
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \todo
    !>  \pmed
    !>  The examples of this generic interface should be extended to corrected weighted covariance matrices.<br>
    !>
    !>  \final{setCovMean}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX

    ! XY - no weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWNO_XY_CK5(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWNO_XY_CK4(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWNO_XY_CK3(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWNO_XY_CK2(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWNO_XY_CK1(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWNO_XY_RK5(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWNO_XY_RK4(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWNO_XY_RK3(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWNO_XY_RK2(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWNO_XY_RK1(cov, mean, x, y, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    ! XY - integer weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWTI_XY_CK5(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWTI_XY_CK4(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWTI_XY_CK3(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWTI_XY_CK2(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWTI_XY_CK1(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWTI_XY_RK5(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWTI_XY_RK4(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWTI_XY_RK3(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWTI_XY_RK2(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWTI_XY_RK1(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    ! XY - real weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWTR_XY_CK5(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWTR_XY_CK4(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWTR_XY_CK3(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWTR_XY_CK2(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWTR_XY_CK1(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: x(:), y(:)
        complex(TKG)            , intent(out)   , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWTR_XY_RK5(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWTR_XY_RK4(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWTR_XY_RK3(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWTR_XY_RK2(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWTR_XY_RK1(cov, mean, x, y, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XY_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: x(:), y(:)
        real(TKG)               , intent(out)   , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! UXD - no weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_CK5(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_CK4(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_CK3(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_CK2(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_CK1(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_RK5(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_RK4(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_RK3(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_RK2(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWNO_UXD_RK1(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    ! UXD - integer weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_CK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_CK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_CK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_CK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_CK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_RK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_RK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_RK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_RK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWTI_UXD_RK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    ! UXD - real weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_CK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_CK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_CK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_CK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_CK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_RK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_RK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_RK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_RK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWTR_UXD_RK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(uppDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! XLD - no weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_CK5(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_CK4(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_CK3(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_CK2(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_CK1(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_RK5(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_RK4(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_RK3(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_RK2(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWNO_XLD_RK1(cov, subset, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWNO_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    ! XLD - integer weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_CK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_CK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_CK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_CK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_CK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_RK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_RK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_RK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_RK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWTI_XLD_RK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTI_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        integer(IK)             , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        integer(IK)             , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

    ! XLD - real weight.

    interface setCovMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_CK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_CK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_CK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_CK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_CK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        complex(TKG)            , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        complex(TKG)            , intent(in)    , contiguous        :: sample(:,:)
        complex(TKG)            , intent(inout) , contiguous        :: cov(:,:)
        complex(TKG)            , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_RK5(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_RK4(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_RK3(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_RK2(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanWTR_XLD_RK1(cov, subset, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanWTR_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                        :: dim
        real(TKG)               , intent(out)                       :: weisum
        real(TKG)               , intent(out)   , contiguous        :: mean(:)
        real(TKG)               , intent(in)    , contiguous        :: weight(:)
        real(TKG)               , intent(in)    , contiguous        :: sample(:,:)
        real(TKG)               , intent(inout) , contiguous        :: cov(:,:)
        real(TKG)               , intent(in)    , contiguous        :: meang(:)
        type(lowDia_type)       , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setCovMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the merged covariance of a sample resulting from the merger of two separate (potentially weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleCov](@ref pm_sampleCov) for more information and definition online updating of sample covariance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** covariance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleCov](@ref pm_sampleCov).<br>
    !>
    !>  \param[in]      covB        :   The input object of the same type and kind and rank and shape as the input argument `covA`,
    !>                                  containing the **biased** covariance of the second sample that must be merged with the first sample.<br>
    !>                                  Only the upper-diagonal triangle of `covB` is referenced.<br>
    !>  \param[in]      covA        :   The input matrix of shape `(1:ndim, 1:ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the **biased** covariance of the first sample that must be merged with the second sample.<br>
    !>                                  Only the upper-diagonal triangle of `covA` is referenced.<br>
    !>  \param[in]      meanDiff    :   The input vector of the same type and kind as the input argument `covA` of shape `(1:ndim)`,
    !>                                  containing the difference of mean of the two samples `meanDiff = meanA - meanB`.<br>
    !>                                  The subtraction order is immaterial.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as `covA`,
    !>                                  containing the sum of the weights of all points in sample \f$A\f$ divided by the sum of the weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>
    !>  \return
    !>  `cov`                       :   The output object of the same type and kind and rank and shape as `meanDiff`,
    !>                                  containing the **full** matrix of **biased** covariance of the sample resulting form the merger of the two samples.<br>
    !>
    !>  \interface{getCovMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: getCovMerged
    !>
    !>      cov(1:ndim, 1:ndim) = getCovMerged(covB(1:ndim, 1:ndim), covA(1:ndim, 1:ndim), meanDiff(1:ndim), fracA)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(covB) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(covA) == size(meanDiff))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{getCovMerged}
    !>  \include{lineno} example/pm_sampleCov/getCovMerged/main.F90
    !>  \compilef{getCovMerged}
    !>  \output{getCovMerged}
    !>  \include{lineno} example/pm_sampleCov/getCovMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \final{getCovMerged}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! RDP_ULD

    interface getCovMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_CK5(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)                                            :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_CK4(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)                                            :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_CK3(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)                                            :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_CK2(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)                                            :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_CK1(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)                                            :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_RK5(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)                                               :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_RK4(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)                                               :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_RK3(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)                                               :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_RK2(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)                                               :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getCovMergedNew_RDP_ULD_RK1(covB, covA, meanDiff, fracA) result(cov)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMergedNew_RDP_ULD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)                                               :: cov(size(covA, 1, IK), size(covA, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the merged covariance of a sample resulting from the merger of two separate (potentially weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleCov](@ref pm_sampleCov) for more information and definition online updating of sample covariance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** covariance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleCov](@ref pm_sampleCov).<br>
    !>
    !>  \param[out]     cov         :   The output object of the same type and kind and rank and shape as `covA`,
    !>                                  containing the **biased** covariance of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged covariance will be written to the argument `covB`.)
    !>  \param[inout]   covB        :   The input or input/output object of the same type and kind and rank and shape as the input argument `covA`,
    !>                                  containing the **biased** covariance of the second sample that must be merged with the first sample.<br>
    !>                                  If the input argument `cov` is missing, then `covB` contains the updated **biased** covariance of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `covB` remain intact upon return.<br>
    !>  \param[in]      covA        :   The input `contiguous` matrix of shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the **biased** covariance of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      meanDiff    :   The input object of the same type and kind as the input argument `covA` of size `size(covA, 1)`,
    !>                                  containing the difference of mean of the two samples `meanDiff = meanA - meanB`.<br>
    !>                                  The subtraction order is immaterial.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as `covA`,
    !>                                  containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>  \param[in]      subset      :   The input scalar constant that can be:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input `cov`, `covB`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input `cov`, `covB`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>
    !>  \interface{setCovMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: setCovMerged
    !>
    !>      call setCovMerged(covB(1:ndim, 1:ndim), covA(1:ndim, 1:ndim), meanDiff(1:ndim), fracA, subset)
    !>      call setCovMerged(cov(1:ndim, 1:ndim), covB(1:ndim, 1:ndim), covA(1:ndim, 1:ndim), meanDiff(1:ndim), fracA, subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(covB) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(meanDiff) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(cov) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{setCovMerged}
    !>  \include{lineno} example/pm_sampleCov/setCovMerged/main.F90
    !>  \compilef{setCovMerged}
    !>  \output{setCovMerged}
    !>  \include{lineno} example/pm_sampleCov/setCovMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! New_RDP_UXD

    interface setCovMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_CK5(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_CK4(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_CK3(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_CK2(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_CK1(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_RK5(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_RK4(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_RK3(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_RK2(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMergedNew_RDP_UXD_RK1(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! New_RDP_XLD

    interface setCovMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_CK5(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_CK4(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_CK3(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_CK2(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_CK1(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_RK5(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_RK4(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_RK3(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_RK2(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMergedNew_RDP_XLD_RK1(cov, covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedNew_RDP_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), covA(:,:), meanDiff(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_RDP_UXD

    interface setCovMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_CK5(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_CK4(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_CK3(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_CK2(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_CK1(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_RK5(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_RK4(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_RK3(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_RK2(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMergedOld_RDP_UXD_RK1(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_RDP_XLD

    interface setCovMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_CK5(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_CK4(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_CK3(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_CK2(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_CK1(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_RK5(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_RK4(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_RK3(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_RK2(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMergedOld_RDP_XLD_RK1(covB, covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMergedOld_RDP_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the merged covariance and mean of a sample resulting from the merger of two separate (potentially weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleCov](@ref pm_sampleCov) for more information and definition online updating of sample covariance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** covariance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleCov](@ref pm_sampleCov).<br>
    !>
    !>  \param[out]     cov         :   The output object of the same type and kind and rank and shape as `covA`,
    !>                                  containing the **biased** covariance of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged covariance will be written to the argument `covB`.)
    !>  \param[out]     mean        :   The output object of the same type and kind and rank and shape as `meanA`,
    !>                                  containing the mean of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged mean will be written to the argument `meanB`.)
    !>  \param[inout]   covB        :   The input or input/output object of the same type and kind and rank and shape as the input argument `covA`,
    !>                                  containing the **biased** covariance of the second sample that must be merged with the first sample.<br>
    !>                                  If the input argument `cov` is missing, then `covB` contains the updated **biased** covariance of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `covB` remain intact upon return.<br>
    !>  \param[in]      meanB       :   The input or input/output object of the same type and kind and rank as the input argument `covA` of size `size(covA, 1)`,
    !>                                  containing the mean of the second sample that must be merged with the mean of the first sample `meanA`.<br>
    !>                                  If the input argument `mean` is missing, then `meanB` contains the updated mean of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `meanB` remain intact upon return.<br>
    !>  \param[in]      covA        :   The input `contiguous` matrix of shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the **biased** covariance of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      meanA       :   The input object of the same type and kind as the input argument `covA` of size `size(covA, 1)`,
    !>                                  containing the mean of the first sample that must be merged with the mean of the second sample `meanB`.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as `covA`,
    !>                                  containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>  \param[in]      subset      :   The input scalar constant that can be:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input `cov`, `covB`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input `cov`, `covB`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>
    !>  \interface{setCovMeanMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: setCovMeanMerged
    !>
    !>      call setCovMeanMerged(covB(1:ndim, 1:ndim), meanB(1:ndim), covA(1:ndim, 1:ndim), meanA(1:ndim), fracA, subset)
    !>      call setCovMeanMerged(cov(1:ndim, 1:ndim), mean(1:ndim), covB(1:ndim, 1:ndim), meanB(1:ndim), covA(1:ndim, 1:ndim), meanA(1:ndim), fracA, subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(covB) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(meanB) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(meanA) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(mean) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(cov) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{setCovMeanMerged}
    !>  \include{lineno} example/pm_sampleCov/setCovMeanMerged/main.F90
    !>  \compilef{setCovMeanMerged}
    !>  \output{setCovMeanMerged}
    !>  \include{lineno} example/pm_sampleCov/setCovMeanMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! New_RDP_UXD

    interface setCovMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_CK5(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_CK4(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_CK3(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_CK2(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_CK1(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_RK5(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_RK4(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_RK3(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_RK2(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_UXD_RK1(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! New_RDP_XLD

    interface setCovMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_CK5(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_CK4(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_CK3(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_CK2(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_CK1(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_RK5(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_RK4(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_RK3(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_RK2(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanMergedNew_RDP_XLD_RK1(cov, mean, covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedNew_RDP_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(out)   , contiguous        :: cov(:,:), mean(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_RDP_UXD

    interface setCovMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_CK5(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_CK4(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_CK3(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_CK2(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_CK1(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_RK5(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_RK4(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_RK3(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_RK2(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_UXD_RK1(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_RDP_XLD

    interface setCovMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_CK5(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_CK4(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_CK3(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_CK2(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_CK1(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        complex(TKG)        , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_RK5(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_RK4(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_RK3(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_RK2(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanMergedOld_RDP_XLD_RK1(covB, meanB, covA, meanA, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanMergedOld_RDP_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(inout) , contiguous        :: covB(:,:), meanB(:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the covAariance resulting from the merger of two separate (potentially weighted) non-singular and singular samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleCov](@ref pm_sampleCov) for more information and definition online updating of sample covAariance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** covAariance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleCov](@ref pm_sampleCov).<br>
    !>
    !>  \param[in]      covA        :   The input/output `contiguous` matrix of shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  
    !>                                  On input, it must contain the **biased** covAariance of the initial non-singular 
    !>                                  sample \f$A\f$ that must be merged with the second singular sample \f$B\f$.<br>
    !>  \param[in]      meanDiff    :   The input object of the same type and kind as the input argument `covA` of size `size(covA, 1)`,
    !>                                  containing the difference of mean of the two non-singular and singular samples `meanDiff = meanA - meanB`.<br>
    !>                                  The subtraction order is immaterial.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as `covA`,
    !>                                  containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>  \param[in]      subset      :   The input scalar constant that can be:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input `covA`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input `covA`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>
    !>  \interface{setCovUpdated}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: setCovUpdated
    !>
    !>      call setCovUpdated(covA(1:ndim, 1:ndim), meanDiff(1:ndim), fracA, subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(meanDiff) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovUpdated](@ref pm_sampleCov::setCovUpdated)<br>
    !>  [setCovMeanMerged](@ref pm_sampleCov::setCovMeanMerged)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{setCovUpdated}
    !>  \include{lineno} example/pm_sampleCov/setCovUpdated/main.F90
    !>  \compilef{setCovUpdated}
    !>  \output{setCovUpdated}
    !>  \include{lineno} example/pm_sampleCov/setCovUpdated/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! Old_RDP_UXD

    interface setCovUpdated

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_CK5(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_CK4(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_CK3(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_CK2(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_CK1(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_RK5(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_RK4(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_RK3(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_RK2(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_UXD_RK1(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_RDP_XLD

    interface setCovUpdated

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_CK5(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_CK4(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_CK3(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_CK2(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_CK1(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:)
        complex(TKG)        , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_RK5(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_RK4(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_RK3(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_RK2(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovUpdatedOld_RDP_XLD_RK1(covA, meanDiff, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovUpdatedOld_RDP_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:)
        real(TKG)           , intent(in)    , contiguous        :: meanDiff(:)
        real(TKG)           , intent(in)                        :: fracA
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the covariance and mean of a sample that results from the merger of two separate (potentially weighted) non-singular \f$A\f$ and singular \f$B\f$ samples.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleCov](@ref pm_sampleCov) for more information and definition online updating of sample covariance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** covariance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleCov](@ref pm_sampleCov).<br>
    !>
    !>  \param[in]      covA        :   The input `contiguous` matrix of shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the **biased** covariance of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      meanA       :   The input object of the same type and kind as the input argument `covA` of size `size(covA, 1)`,
    !>                                  containing the mean of the first sample that must be merged with the mean of the second sample `meanB`.<br>
    !>  \param[in]      meanB       :   The input or input/output object of the same type and kind and rank as the input argument `covA` of size `size(covA, 1)`,
    !>                                  containing the mean of the second sample that must be merged with the mean of the first sample `meanA`.<br>
    !>                                  If the input argument `mean` is missing, then `meanB` contains the updated mean of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `meanB` remain intact upon return.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as `covA`,
    !>                                  containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>  \param[in]      subset      :   The input scalar constant that can be:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [lowDia](@ref pm_matrixSubset::lowDia), implying that only the lower-diagonal subset of the input `cov`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                      <li>    The constant [uppDia](@ref pm_matrixSubset::uppDia), implying that only the upper-diagonal subset of the input `cov`, and `covA` matrices must be accessed and/or manipulated.<br>
    !>                                  </ol>
    !>                                  This input argument is merely serves to resolve the different procedures of this generic interface from each other at compile-time.<br>
    !>
    !>  \interface{setCovMeanUpdated}
    !>  \code{.F90}
    !>
    !>      use pm_sampleCov, only: setCovMeanUpdated
    !>
    !>      call setCovMeanUpdated(covA(1:ndim, 1:ndim), meanA(1:ndim), meanB(1:ndim), fracA, subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(meanB) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(meanA) == shape(covA))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [setCovMeanUpdated](@ref pm_sampleCov::setCovMeanUpdated)<br>
    !>  [getVarCorrection](@ref pm_sampleVar::getVarCorrection)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>
    !>  \example{setCovMeanUpdated}
    !>  \include{lineno} example/pm_sampleCov/setCovMeanUpdated/main.F90
    !>  \compilef{setCovMeanUpdated}
    !>  \output{setCovMeanUpdated}
    !>  \include{lineno} example/pm_sampleCov/setCovMeanUpdated/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleCov](@ref test_pm_sampleCov)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

    ! New_RDP_UXD

    interface setCovMeanUpdated

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_CK5(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_CK4(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_CK3(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_CK2(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_CK1(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_RK5(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_RK4(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_RK3(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_RK2(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_UXD_RK1(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_UXD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(uppDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! New_RDP_XLD

    interface setCovMeanUpdated

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_CK5(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_CK4(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_CK3(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_CK2(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_CK1(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)           , intent(in)                        :: fracA
        complex(TKG)        , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        complex(TKG)        , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_RK5(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_RK4(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_RK3(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_RK2(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCovMeanUpdatedOld_RDP_XLD_RK1(covA, meanA, meanB, fracA, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCovMeanUpdatedOld_RDP_XLD_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)                        :: fracA
        real(TKG)           , intent(inout) , contiguous        :: covA(:,:), meanA(:)
        real(TKG)           , intent(in)    , contiguous        :: meanB(:)
        type(lowDia_type)   , intent(in)                        :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleCov ! LCOV_EXCL_LINE