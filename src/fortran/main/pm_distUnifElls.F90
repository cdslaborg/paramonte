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
!>  This module contains classes and procedures for computing various statistical
!>  quantities related to the <b>Multiple MultiVariate Uniform Ellipsoid (MMVUE) distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following
!>  quantities of the <b>Multiple MultiVariate Uniform Ellipsoid (MMVUE) distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>  </ol>
!>
!>  An \f$\ndim\f$-dimensional MMVUE distribution is represented by a set of \f$\ndim\f$-dimensional hyper-ellipsoid supports.<br>
!>  The MMVUE distribution is fully determined by a set of ellipsoids \f$\ell_i, i = 1 : N_\ell\f$ containing its support.<br>
!>  An ellipsoid \f$\ell\f$ is in turn fully determined by its representative Gramian matrix \f$\gramian_\ell\f$,
!>  containing all points in \f$\mathbb{R}^\ndim\f$ that satisfy,
!>  \f{equation}{
!>      \large
!>      (X - \mu_\ell)^T ~ \gramian_\ell^{-1} ~ (X - \mu_\ell) \leq 1 ~,
!>  \f}
!>  where \f$\mu_\ell\f$ represents the center of the ellipsoid,
!>  \f$(X - \mu_\ell)^T\f$ is the transpose of the vector \f$(X - \mu_\ell)\f$,
!>  and \f$\gramian_\ell^{-1}\f$ is the inverse of the matrix \f$\gramian_\ell\f$.<br>
!>
!>  The volume of this ellipsoid is given by,
!>  \f{equation}{
!>      \large
!>      V(\ell) = V_\ndim \sqrt{\left| \gramian_\ell \right|} ~,
!>  \f}
!>  where \f$\left|\gramian_\ell\right|\f$ is the determinant of \f$\gramian_\ell\f$ and,
!>  \f{equation}{
!>      \large
!>      V_\ndim =
!>      \frac{\pi^{\ndim / 2}}{\up\Gamma(1 + \ndim / 2)} =
!>      \begin{cases}
!>          \frac{1}{(\ndim/2)!} \pi^{\ndim/2} & \text{if $\ndim$ is even} \\\\
!>          2^\ndim \frac{1}{\ndim!} \big( \frac{\ndim-1}{2} \big)! ~ \pi^{(\ndim-1)/2} & \text{if $\ndim$ is odd}
!>      \end{cases}
!>  \f}
!>  is the volume of an \f$\ndim\f$-ball (that is, a unit-logChoDia \f$\ndim\f$-dimensional hyper-sphere).<br>
!>  It is readily seen that the corresponding unit-volume ellipsoid \f$\widehat\ell\f$ has the representative Gramian matrix,
!>  \f{equation}{
!>      \large
!>      \gramian_{\widehat\ell} = V_\ell^{-2/\ndim} ~ \gramian_\ell ~.
!>  \f}
!>  More generally, to scale the ellipsoidal support of an MMVUE distribution \f$\ell\f$ by some factor \f$\alpha\f$
!>  along each coordinate axis, it suffices to be used the new scaled ellipsoid \f$\ell^*\f$ with the representative Gramian matrix,
!>  \f{equation}{
!>      \large
!>      \gramian_{\ell^*} = \alpha^2 ~ \gramian_\ell ~,
!>  \f}
!>  in which case, the volume of \f$\ell^*\f$ becomes,
!>  \f{equation}{
!>      \large
!>      V_{\ell^*} = \alpha^\ndim ~ V_\ell ~.
!>  \f}
!>
!>  The Probability Density Function (**PDF**) of the MMVUE distribution
!>  with multiple-ellipsoidal support \f$\ell\f$ is given by,<br>
!>  \f{equation}{
!>      \large
!>      \pi(X | \ell) = \frac{1}{V(\ell)} = \frac{1}{V_\ndim \sqrt{\left|\gramian_\ell\right|}} ~,
!>  \f}
!>  where \f$\sqrt{\left|\gramian_\ell\right|}\f$ represents the union of the volumes of all ellipsoids representing the distribution.<br>
!>  Because there is no closed-form expression for the shared volume of two more overlapping ellipsoids, this union of volumes must be approximated.<br>
!>  Within this module, this approximation is done using a Monte Carlo approach via the generic interface [getUnifEllsLogPDF](@ref pm_distUnifElls::getUnifEllsLogPDF).<br>
!>
!>  **Random Number Generation**<br>
!>
!>  Note that in an MMVUE distribution, multiple ellipsoids can have partially overlapping support,
!>  in which case, the shared support is counted only once in random number generation.<br>
!>  Therefore, the random number generation has to be carefully done so that overlapping ellipsoids
!>  are not sampled too often.<br>
!>
!>  The **RNG** generic interfaces within this module generate uniformly-distributed random vectors from within a single \f$\ndim\f$-dimensional
!>  hyper-ellipsoid by generalizing the proposed approach of Marsaglia (1972) for choosing a point from the surface of a sphere.<br>
!>  <ol>
!>      <li>    Generate a normalized (unit) \f$\ndim\f$-dimensional [Multivariate Normal random vector](@ref pm_distMultiNorm),
!>              \f{equation}{
!>                  \unit{r} = \frac{1}{\sqrt{\sum_{i=1}^\ndim g_i^2}} \sum_{j=1}^{\ndim} g_j \unit{r}_j ~,
!>              \f}
!>              where \f$\unit{r}_j\f$ are the unit vectors representing the orthogonal basis the of \f$\ndim\f$-space.<br>
!>              This unit vector \f$\unit{r}\f$ has uniform random orientation and distribution on the surface of an \f$\ndim\f$-dimensional
!>              sphere of unit radius centered at the origin of the coordinate system.<br>
!>      <li>    Generate a uniformly-distributed random number \f$u\in\mathcal{U}[0,1)\f$ so that,
!>              \f{equation}{
!>                  \bs{r}_{\sphere} = u^{1/\ndim} \unit{r} ~,
!>              \f}
!>              represents a vector pointing to a uniformly-distributed location inside of the \f$\ndim\f$-sphere.<br>
!>      <li>    To obtain a uniformly-distributed vector from inside of an \f$\ndim\f$-dimensional ellipsoid,
!>              first compute the [Cholesky decomposition](@ref pm_matrixChol) of the representative Gramian matrix \f$\gramian\f$ of the ellipsoid,
!>              \f{equation}{
!>                  \gramian = \mat{L}\mat{L}^H ~,
!>              \f}
!>              where \f$\mat{L}\f$ is the left triangular matrix resulting from the Cholesky factorization,
!>              and \f$\mat{L}^T\f$ is its Hermitian transpose.<br>
!>              Then the vector,
!>              \f{equation}{
!>                  \bs{r}_\ell = \mat{L} ~ \bs{r}_\sphere + \mu_\ell ~,
!>              \f}
!>              is uniformly distributed inside ellipsoid \f$\ell\f$ centered at \f$\mu_\ell\f$.<br>
!>  </ol>
!>
!>  \note
!>  <ol>
!>      <li>    An ellipsoid with a diagonal representative Gramian matrix \f$\gramian_\ell\f$ is an
!>              uncorrelated **hyper-ellipsoid** whose axes are parallel to the coordinate axes.<br>
!>      <li>    An uncorrelated hyper-ellipsoid with equal diagonal elements in its
!>              representative Gramian matrix \f$\gramian_\ell\f$ is a **hyper-sphere**.<br>
!>  </ol>
!>
!>  <b>The covariance matrix of a single ellipsoidal support in the MMVUE distribution</b><br>
!>  The covariance matrix of the multivariate uniform ellipsoidal distribution is given by its Gramian matrix as,
!>  \f{eqnarray*}{
!>       \Sigma = \frac{\gramian}{\ndim + 2} ~,
!>  \f}
!>
!>  \see
!>  [pm_ellipsoid](@ref pm_ellipsoid)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>  [pm_distUnifElls](@ref pm_distUnifElls)<br>
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  Marsaglia G., et al. (1972), Choosing a point from the surface of a sphere. The Annals of Mathematical Statistics 43(2):645-646.<br>
!>  Gammell JD, Barfoot TD (2014) The probability density function of a transformation-based hyper-ellipsoid sampling technique.<br>
!>
!>  \test
!>  [test_pm_distUnifElls](@ref test_pm_distUnifElls)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distUnifElls

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distUnifElls"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Multiple MultiVariate Uniform Ellipsoid (MMVUE)
    !>  as defined in the description of [pm_distUnifElls](@ref pm_distUnifElls).
    !>
    !>  \details
    !>  See the documentation of [pm_distUnifElls](@ref pm_distUnifElls) for the definition of the Multiple MultiVariate Uniform Ellipsoid (MMVUE) distribution.
    !>
    !>  \interface{distUnifElls_type}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifElls, only: distUnifElls_type
    !>      type(distUnifElls_type) :: distUnifEll
    !>
    !>      distUnifEll = distUnifElls_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distUnifElls](@ref test_pm_distUnifElls)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distUnifElls_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distUnifElls_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of an approximation of the Probability Density Function (PDF)
    !>  of the Multiple MultiVariate Uniformly Ellipsoidal (MMVUE) Distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distUnifElls](@ref pm_distUnifElls) for details of the definition of the PDF.<br>
    !>
    !>  \note
    !>  Note that the procedures of this generic interface do not require an
    !>  input `X` value representing a location within the domain of the density function.<br>
    !>  This is fine and intentional by design because the PDF is uniform across the entire support of the PDF.<br>
    !>
    !>  \param[in]      invmul      :   The input vector of shape `(1:nsam)` of the same type and kind as the output `logPDF`,
    !>                                  each element of which represents the inverse of the number of ellipsoids within which the corresponding sampled random vector falls.
    !>  \param[in]      ndim        :   The input positive scalar of type `integer` of default kind \IK, containing the number of dimensions of the domain of the distribution.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `nell` is also present and `chol` is missing.)
    !>  \param[in]      nell        :   The input positive scalar of type `integer` of default kind \IK, containing the number of ellpisoids in the distribution.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `ndim` is also present and `chol` is missing.)
    !>  \param[in]      chol        :   The input `contiguous` array of shape `(1:ndim, 1:ndim, 1:nell)` whose specified triangular `subset` contains
    !>                                  the [Cholesky Factorization](@ref pm_matrixChol) of the Gramian matrix of the MMVUE distribution.<br>
    !>                                  (**optional**, the default is the Identity matrix of rank `ndim`. It must be present <b>if and only if</b> the input argument `subset` and `invGram` are also present.)
    !>
    !>  \return
    !>  `logPDF`                    :   The output scalar of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the natural logarithm of the PDF of the distribution.
    !>
    !>  \interface{getUnifEllsLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifElls, only: getUnifEllsLogPDF
    !>
    !>      logPDF = getUnifEllsLogPDF(chol, invmul)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < nell` must hold for the corresponding input arguments.<br>
    !>  The condition `all(invmul <= 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(chol, 1) == size(chol, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when the input argument `ndim` is present.<br>
    !>
    !>  \see
    !>  [pm_ellipsoid](@ref pm_ellipsoid)<br>
    !>  [getUnifEllsRand](@ref pm_distUnifElls::getUnifEllsRand)<br>
    !>  [setUnifEllsRand](@ref pm_distUnifElls::setUnifEllsRand)<br>
    !>
    !>  \example{getUnifEllsLogPDF}
    !>  \include{lineno} example/pm_distUnifElls/getUnifEllsLogPDF/main.F90
    !>  \compilef{getUnifEllsLogPDF}
    !>  \output{getUnifEllsLogPDF}
    !>  \include{lineno} example/pm_distUnifElls/getUnifEllsLogPDF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnifElls](@ref test_pm_distUnifElls)
    !>
    !>  \final{getUnifEllsLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! NELL

    interface getUnifEllsLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUnifEllsLogPDFNELL_RK5(invmul, ndim, nell) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFNELL_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(in)                    :: ndim, nell
        real(RKG)           , intent(in)    , contiguous    :: invmul(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getUnifEllsLogPDFNELL_RK4(invmul, ndim, nell) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFNELL_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(in)                    :: ndim, nell
        real(RKG)           , intent(in)    , contiguous    :: invmul(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getUnifEllsLogPDFNELL_RK3(invmul, ndim, nell) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFNELL_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(in)                    :: ndim, nell
        real(RKG)           , intent(in)    , contiguous    :: invmul(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getUnifEllsLogPDFNELL_RK2(invmul, ndim, nell) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFNELL_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(in)                    :: ndim, nell
        real(RKG)           , intent(in)    , contiguous    :: invmul(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getUnifEllsLogPDFNELL_RK1(invmul, ndim, nell) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFNELL_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(in)                    :: ndim, nell
        real(RKG)           , intent(in)    , contiguous    :: invmul(:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CHOL

    interface getUnifEllsLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUnifEllsLogPDFCHOL_RK5(invmul, chol) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFCHOL_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: invmul(:), chol(:,:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getUnifEllsLogPDFCHOL_RK4(invmul, chol) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFCHOL_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: invmul(:), chol(:,:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getUnifEllsLogPDFCHOL_RK3(invmul, chol) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFCHOL_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: invmul(:), chol(:,:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getUnifEllsLogPDFCHOL_RK2(invmul, chol) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFCHOL_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: invmul(:), chol(:,:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getUnifEllsLogPDFCHOL_RK1(invmul, chol) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllsLogPDFCHOL_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: invmul(:), chol(:,:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a collection of random vectors of size `ndim` from the `ndim`-dimensional Multiple MultiVariate Uniform Ellipsoidal (MMVUE) distribution,
    !>  with the specified input `mean(1:ndim, 1:nell)` and optionally the specified `subset` of the Cholesky Factorization of the Gramian matrices of the MMVUE distribution.<br>
    !>
    !>  \details
    !>  Here,
    !>  <ol>
    !>      <li>    the variable `ndim` represents the number of dimensions of the domain of the distribution (that is, the number of dimensions of the ellipsoids),
    !>      <li>    the variable `nsam` represents the number of randomly sampled points from the set of ellipsoids specified by the input arguments,
    !>      <li>    the variable `nell` represents the number of ellipsoids in the distirbution.
    !>  </ol>
    !>
    !>  \param[inout]   rng         :   The input/output scalar that can be an object of,
    !>                                  <ol>
    !>                                      <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                              implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                      <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                              implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                                  </ol>
    !>  \param[out]     rand        :   The output `contiguous` matrix of shape `(1:ndim, 1:nsam)` of<br>
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the random output vectors.<br>
    !>  \param[out]     mahalSq     :   The output matrix of shape `(1:nell, 1:nsam)` of the same type and kind as the output `rand`,
    !>                                  containing the squared Mahalanobis distances of individual randomly sampled vectors `rand(:, isam)` from the
    !>                                  centers of the specified ellipsoids in the distribution via the input argument `mean(1:ndim, 1:nell)`.<br>
    !>                                  A value larger than `1` for `mahalSq(iell, isam)` implies that the ranodm vector `rand(1:ndim, isam)` is
    !>                                  outside ellipsoid `iell` specified by the input arguments, while a value less than `1` implies that the
    !>                                  random vector is within the ellipsoid.<br>
    !>                                  Thus, the expression `count(mahalSq(1:nell, isam) <= 1)` yields the membership count of `rand(1:ndim, isam)` in all input ellipsoids.
    !>  \param[out]     invmul      :   The output vector of shape `(1:nsam)` of the same type and kind as the output `rand`,
    !>                                  each element of which represents the inverse of the number of ellipsoids within which the corresponding random vector `rand(1:ndim, isam)` falls.
    !>  \param[in]      mean        :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`,
    !>                                  representing the centers of the ellipsoids in the distribution.<br>
    !>  \param[in]      chol        :   The input `contiguous` array of shape `(1:ndim, 1:ndim, 1:nell)` whose specified triangular `subset` contains
    !>                                  the [Cholesky Factorization](@ref pm_matrixChol) of the Gramian matrix of the MMVUE distribution.<br>
    !>                                  (**optional**, the default is the Identity matrix of rank `ndim`. It must be present <b>if and only if</b> the input argument `subset` and `invGram` are also present.)
    !>  \param[in]      subset      :   The input scalar constant that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                              implying that the upper-diagonal triangular block of the input `chol` must be used while the lower subset is not referenced.<br>
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                              implying that the lower-diagonal triangular block of the input `chol` must be used while the upper subset is not referenced.<br>
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `chol` is present.)
    !>  \param[in]      invGram     :   The input array of shape `(1:ndim, 1:ndim, 1:nell)` of the same type and kind as the input argument `point`,
    !>                                  containing the collection of square matrices of full inverse representative Gramian matrix of the \f$\ndim\f$-dimensional ellipsoid.<br>
    !>                                  This argument is needed to determine the membership of the output random vectors in the ellipsoids.<br>
    !>                                  (**optional**. It must be present **if and only if** the inpt argument `chol` is present.)
    !>
    !>  \interface{setUnifEllsRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifElls, only: setUnifEllsRand
    !>
    !>      call setUnifEllsRand(rng, rand(1:ndim, 1:nsam), mahalSq(1:nell, 1:nsam), invmul(1:nsam), mean(1:ndim, 1:nell))
    !>      call setUnifEllsRand(rng, rand(1:ndim, 1:nsam), mahalSq(1:nell, 1:nsam), invmul(1:nsam), mean(1:ndim, 1:nell), chol(1:ndim, 1:ndim, 1:nell), subset, invGram(1:ndim, 1:ndim, 1:nell))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mean, 1) == size(rand, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(mahalSq) == [size(mean, 2), size(rand, 2)])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(chol) == [size(mean, 1), size(mean, 1), size(mean, 2)])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(chol) == shape(invGram))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>  The procedures of this generic interface are `pure` when the input argument `rng` is set to
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) and the compile-time macro `CHECK_ENABLED` is set to `0` or is undefined.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [pm_distNorm](@ref pm_distNorm)<br>
    !>  [pm_distUnif](@ref pm_distUnif)<br>
    !>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
    !>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
    !>  [pm_distUnifElls](@ref pm_distUnifElls)<br>
    !>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
    !>
    !>  \example{setUnifEllsRand}
    !>  \include{lineno} example/pm_distUnifElls/setUnifEllsRand/main.F90
    !>  \compilef{setUnifEllsRand}
    !>  \output{setUnifEllsRand}
    !>  \include{lineno} example/pm_distUnifElls/setUnifEllsRand/main.out.F90
    !>  \postproc{setUnifEllsRand}
    !>  \include{lineno} example/pm_distUnifElls/setUnifEllsRand/main.py
    !>  \vis{setUnifEllsRand}
    !>  \image html pm_distUnifElls/setUnifEllsRand/setUnifEllsRandMean.png width=700
    !>  \image html pm_distUnifElls/setUnifEllsRand/setUnifEllsRandMeanChol.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifElls](@ref test_pm_distUnifElls)
    !>
    !>  \final{setUnifEllsRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D2 RNGF

    interface setUnifEllsRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMMUR_RNGF_AM_DC_XXX_D2_RK5(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMMUR_RNGF_AM_DC_XXX_D2_RK4(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMMUR_RNGF_AM_DC_XXX_D2_RK3(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMMUR_RNGF_AM_DC_XXX_D2_RK2(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMMUR_RNGF_AM_DC_XXX_D2_RK1(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_UXD_D2_RK5(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_UXD_D2_RK4(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_UXD_D2_RK3(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_UXD_D2_RK2(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_UXD_D2_RK1(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_XLD_D2_RK5(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_XLD_D2_RK4(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_XLD_D2_RK3(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_XLD_D2_RK2(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMMUR_RNGF_AM_AC_XLD_D2_RK1(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGF_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGX

    interface setUnifEllsRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_DC_XXX_D2_RK5(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_DC_XXX_D2_RK4(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_DC_XXX_D2_RK3(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_DC_XXX_D2_RK2(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_DC_XXX_D2_RK1(rng, rand, mahalSq, invmul, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_UXD_D2_RK5(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_UXD_D2_RK4(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_UXD_D2_RK3(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_UXD_D2_RK2(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_UXD_D2_RK1(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_XLD_D2_RK5(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_XLD_D2_RK4(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_XLD_D2_RK3(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_XLD_D2_RK2(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMMUR_RNGX_AM_AC_XLD_D2_RK1(rng, rand, mahalSq, invmul, mean, chol, subset, invGram)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMMUR_RNGX_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(out)   , contiguous    :: mahalSq(:,:)
        real(RKG)               , intent(out)   , contiguous    :: invmul(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:,:), chol(:,:,:), invGram(:,:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distUnifElls