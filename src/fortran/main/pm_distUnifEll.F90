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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>MultiVariate Uniform Ellipsoid (MVUE) distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>MultiVariate Uniform Ellipsoid (MVUE) distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>  </ol>
!>
!>  An \f$\ndim\f$-dimensional MVUE distribution is represented by an \f$\ndim\f$-dimensional hyper-ellipsoid support.<br>
!>  The MVUE distribution is fully determined by the ellipsoid \f$\ell\f$ containing its support.<br>
!>  The ellipsoid \f$\ell\f$ is in turn fully determined by its representative Gramian matrix \f$\gramian_\ell\f$,
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
!>  More generally, to scale the ellipsoidal support of an MVUE distribution \f$\ell\f$ by some factor \f$\alpha\f$
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
!>  The Probability Density Function (**PDF**) of the MVUE distribution
!>  with ellipsoidal support \f$\ell\f$ is given by,<br>
!>  \f{equation}{
!>      \large
!>      \pi(X | \ell) = \frac{1}{V(\ell)} = \frac{1}{V_\ndim \sqrt{\left|\gramian_\ell\right|}} ~.
!>  \f}
!>
!>  **Random Number Generation**<br>
!>
!>  The **RNG** generic interfaces within this module generate uniformly-distributed random vectors from within an \f$\ndim\f$-dimensional
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
!>  <b>The covariance matrix of the MVUE distribution</b><br>
!>  The covariance matrix of the multivariate uniform ellipsoidal distribution is given by its Gramian matrix as,
!>  \f{eqnarray*}{
!>       \Sigma = \frac{\gramian}{\ndim + 2} ~,
!>  \f}
!>
!>  \see
!>  [pm_ellipsoid](@ref pm_ellipsoid)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>  Marsaglia G., et al. (1972), Choosing a point from the surface of a sphere. The Annals of Mathematical Statistics 43(2):645-646.<br>
!>  Gammell JD, Barfoot TD (2014) The probability density function of a transformation-based hyper-ellipsoid sampling technique.<br>
!>
!>  \test
!>  [test_pm_distUnifEll](@ref test_pm_distUnifEll)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distUnifEll

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distUnifEll"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type MultiVariate Uniform Ellipsoid (MVUE)
    !>  as defined in the description of [pm_distUnifEll](@ref pm_distUnifEll).
    !>
    !>  \details
    !>  See the documentation of [pm_distUnifEll](@ref pm_distUnifEll) for the definition of the MultiVariate Uniform Ellipsoid (MVUE) distribution.
    !>
    !>  \interface{distUnifEll_type}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifEll, only: distUnifEll_type
    !>      type(distUnifEll_type) :: distUnifEll
    !>
    !>      distUnifEll = distUnifEll_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distUnifEll](@ref test_pm_distUnifEll)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distUnifEll_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distUnifEll_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF)
    !>  of the MultiVariate MVUE (MVUE) Distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distUnifEll](@ref pm_distUnifEll) for details of the definition of the PDF.<br>
    !>
    !>  \note
    !>  Note that the procedures of this generic interface do not require an
    !>  input `X` value representing a location within the domain of the density function.<br>
    !>  This is fine and intentional by design because the PDF is uniform across the entire support of the PDF.<br>
    !>
    !>  \param[in]  logChoDia   :   The input scalar or `contiguous` vector of shape `(1:ndim)` of the same type and kind as the output `logPDF`.<br>
    !>                              <ol>
    !>                                  <li>    If `logChoDia` is scalar, it must contain the natural logarithm of the radius
    !>                                          of the corresponding hyper-spherical support of the distribution.<br>
    !>                                  <li>    If `logChoDia` is vector, it must contain the natural logarithm of the diagonal elements of
    !>                                          the Cholesky Factorization of the symmetric positive-definite representative Gramian matrix \f$\gramian_\ell\f$
    !>                                          of the corresponding hyper-spherical support of the distribution.<br>
    !>                              </ol>
    !>                              (**optional**. It must be present **if and only if** the input `gramian` is missing.)
    !>  \param[in]  ndim        :   The input positive scalar of type `integer` of default kind \IK, containing the number of dimensions of the domain of the distribution.<br>
    !>                              (**optional**. It must be present **if and only if** the input `logChoDia` is present and is a scalar.)
    !>  \param[in]  gramian     :   The input square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the output `logPDF`,
    !>                              containing the **upper** triangle and diagonal elements of the
    !>                              [representative Gramian matrix](@ref pm_distUnifEll) \f$\gramian_\ell\f$ of the ellipsoid.<br>
    !>                              (**optional**. It must be present **if and only if** the input `logChoDia` is missing.)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the natural logarithm of the PDF of the distribution.
    !>
    !>  \interface{getUnifEllLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifEll, only: getUnifEllLogPDF
    !>
    !>      logPDF = getUnifEllLogPDF(logChoDia, ndim)         ! Hyper-sphere support
    !>      logPDF = getUnifEllLogPDF(logChoDia(1:ndim))       ! Hyper-ellipsoid support
    !>      logPDF = getUnifEllLogPDF(gramian(1:ndim, 1:ndim)) ! Hyper-ellipsoid support
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold for the corresponding input arguments.<br>
    !>  The condition `size(gramian, 1) == size(gramian, 2)` must hold for the corresponding input arguments.<br>
    !>  The condition \f$0 < |\ms{gramian}|\f$ must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when the input argument `ndim` is present.<br>
    !>
    !>  \see
    !>  [pm_ellipsoid](@ref pm_ellipsoid)<br>
    !>  [getUnifEllRand](@ref pm_distUnifEll::getUnifEllRand)<br>
    !>  [setUnifEllRand](@ref pm_distUnifEll::setUnifEllRand)<br>
    !>
    !>  \example{getUnifEllLogPDF}
    !>  \include{lineno} example/pm_distUnifEll/getUnifEllLogPDF/main.F90
    !>  \compilef{getUnifEllLogPDF}
    !>  \output{getUnifEllLogPDF}
    !>  \include{lineno} example/pm_distUnifEll/getUnifEllLogPDF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distUnifEll](@ref test_pm_distUnifEll)
    !>
    !>  \final{getUnifEllLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getUnifEllLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getUnifEllLogPDF_D0_RK5(logChoDia, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logChoDia
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getUnifEllLogPDF_D0_RK4(logChoDia, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logChoDia
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getUnifEllLogPDF_D0_RK3(logChoDia, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logChoDia
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getUnifEllLogPDF_D0_RK2(logChoDia, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logChoDia
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getUnifEllLogPDF_D0_RK1(logChoDia, ndim) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                            :: ndim
        real(RKG)   , intent(in)                            :: logChoDia
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUnifEllLogPDF_D1_RK5(logChoDia) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: logChoDia(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getUnifEllLogPDF_D1_RK4(logChoDia) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: logChoDia(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getUnifEllLogPDF_D1_RK3(logChoDia) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: logChoDia(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getUnifEllLogPDF_D1_RK2(logChoDia) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: logChoDia(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getUnifEllLogPDF_D1_RK1(logChoDia) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: logChoDia(:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getUnifEllLogPDF_D2_RK5(gramian) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: gramian(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getUnifEllLogPDF_D2_RK4(gramian) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: gramian(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getUnifEllLogPDF_D2_RK3(gramian) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: gramian(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getUnifEllLogPDF_D2_RK2(gramian) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: gramian(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getUnifEllLogPDF_D2_RK1(gramian) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUnifEllLogPDF_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: gramian(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a (collection) of random vector(s) of size `ndim` from the `ndim`-dimensional MultiVariate Uniform Ellipsoidal (MVUE) distribution,
    !>  optionally with the specified input `mean(1:ndim)` and the specified `subset` of the Cholesky Factorization of the Gramian matrix of the MVUE distribution.
    !>
    !>  \details
    !>  The procedures of this generic interface are merely wrappers around the subroutine interface [setUnifRand](@ref pm_distUnifEll::setUnifEllRand).<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type).)
    !>  \param[in]      mean    :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`, representing the mean of the MVUE distribution.<br>
    !>                              (**optional**, default = `[(0., i = 1, size(rand))]`. It must be present if the input argument `chol` is missing.)
    !>  \param[in]      chol    :   The input `contiguous` matrix of shape `(ndim, ndim)` whose specified triangular `subset` contains the [Cholesky Factorization](@ref pm_matrixChol) of the Gramian matrix of the MVUE distribution.<br>
    !>                              (**optional**, the default is the Identity matrix of rank `ndim`. It must be present <b>if and only if</b> the input argument `subset` is also present.)
    !>  \param[in]      subset  :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                          implying that the upper-diagonal triangular block of the input `chol` must be used while the lower subset is not referenced.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                          implying that the lower-diagonal triangular block of the input `chol` must be used while the upper subset is not referenced.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `chol` is present.)
    !>  \param[in]      nsam    :   The input scalar `integer` of default kind \IK containing the number of random MVUE vectors to generate.<br>
    !>                              (**optional**. If present, the output `rand` is of rank `2`, otherwise is of rank `1`.)
    !>
    !>  \return
    !>  `rand`                  :   The output vector of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the MVUE distributed random output vector(s):<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `nsam` is missing, then `rand` shall be of shape `(1:ndim)`.<br>
    !>                                  <li>    If the input argument `nsam` is present, then `rand` shall be of shape `(1:ndim, 1:nsam)`.<br>
    !>                              </ol>
    !>
    !>  \interface{getUnifEllRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifEll, only: getUnifEllRand
    !>
    !>      ! single vector, using default rng
    !>
    !>      rand(1:ndim) = getUnifEllRand(mean(1:ndim))
    !>      rand(1:ndim) = getUnifEllRand(chol(1:ndim, 1:ndim), subset)
    !>      rand(1:ndim) = getUnifEllRand(mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! single vector, using custom rng
    !>
    !>      rand(1:ndim) = getUnifEllRand(rng, mean(1:ndim))
    !>      rand(1:ndim) = getUnifEllRand(rng, chol(1:ndim, 1:ndim), subset)
    !>      rand(1:ndim) = getUnifEllRand(rng, mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using default rng
    !>
    !>      rand(1:ndim, 1:nsam) = getUnifEllRand(mean(1:ndim), nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifEllRand(chol(1:ndim, 1:ndim), subset, nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifEllRand(mean(1:ndim), chol(1:ndim, 1:ndim), subset, nsam)
    !>
    !>      ! collection of `nsam` vectors, using custom rng
    !>
    !>      rand(1:ndim, 1:nsam) = getUnifEllRand(rng, mean(1:ndim), nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifEllRand(rng, chol(1:ndim, 1:ndim), subset, nsam)
    !>      rand(1:ndim, 1:nsam) = getUnifEllRand(rng, mean(1:ndim), chol(1:ndim, 1:ndim), subset, nsam)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [pm_ellipsoid](@ref pm_ellipsoid)<br>
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>
    !>  \example{getUnifEllRand}
    !>  \include{lineno} example/pm_distUnifEll/getUnifEllRand/main.F90
    !>  \compilef{getUnifEllRand}
    !>  \output{getUnifEllRand}
    !>  \include{lineno} example/pm_distUnifEll/getUnifEllRand/main.out.F90
    !>  \postproc{getUnifEllRand}
    !>  \include{lineno} example/pm_distUnifEll/getUnifEllRand/main.py
    !>  \vis{getUnifEllRand}
    !>  \image html pm_distUnifEll/getUnifEllRand/getUnifEllRandMean.RK.png width=700
    !>  \image html pm_distUnifEll/getUnifEllRand/getUnifEllRandChol.RK.png width=700
    !>  \image html pm_distUnifEll/getUnifEllRand/getUnifEllRandMeanChol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifEll](@ref test_pm_distUnifEll)
    !>
    !>  \final{getUnifEllRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D1 RNGD

    interface getUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D1_RK5(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D1_RK4(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D1_RK3(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D1_RK2(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D1_RK1(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D1_RK5(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D1_RK4(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D1_RK3(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D1_RK2(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D1_RK1(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D1_RK5(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D1_RK4(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D1_RK3(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D1_RK2(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D1_RK1(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D1_RK5(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D1_RK4(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D1_RK3(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D1_RK2(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D1_RK1(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D1_RK5(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D1_RK4(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D1_RK3(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D1_RK2(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D1_RK1(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGF

    interface getUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D1_RK5(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D1_RK4(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D1_RK3(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D1_RK2(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D1_RK1(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGX

    interface getUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D1_RK5(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D1_RK4(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D1_RK3(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D1_RK2(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D1_RK1(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGD

    interface getUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D2_RK5(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D2_RK4(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D2_RK3(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D2_RK2(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_AM_DC_XXX_D2_RK1(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D2_RK5(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D2_RK4(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D2_RK3(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D2_RK2(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_DM_AC_UXD_D2_RK1(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D2_RK5(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D2_RK4(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D2_RK3(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D2_RK2(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_DM_AC_XLD_D2_RK1(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D2_RK5(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D2_RK4(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D2_RK3(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D2_RK2(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_AM_AC_UXD_D2_RK1(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D2_RK5(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D2_RK4(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D2_RK3(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D2_RK2(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGD_AM_AC_XLD_D2_RK1(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGD_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGF

    interface getUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D2_RK5(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D2_RK4(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D2_RK3(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D2_RK2(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_AM_DC_XXX_D2_RK1(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_DM_AC_UXD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_DM_AC_XLD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_AM_AC_UXD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGF_AM_AC_XLD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGF_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGX

    interface getUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D2_RK5(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D2_RK4(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D2_RK3(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D2_RK2(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_AM_DC_XXX_D2_RK1(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_DM_AC_UXD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_DM_AC_XLD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_AM_AC_UXD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMUR_RNGX_AM_AC_XLD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMUR_RNGX_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a (collection) of random vector(s) of size `ndim` from the `ndim`-dimensional MultiVariate Uniform Ellipsoidal (MVUE) distribution,
    !>  optionally with the specified input `mean(1:ndim)` and the specified `subset` of the Cholesky Factorization of the Gramian matrix of the MVUE distribution.
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type).)
    !>  \param[out]     rand    :   The output `contiguous` vector of shape `(1:ndim)` or matrix of shape `(1:ndim, 1:nsam)` of<br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the random output vector(s).<br>
    !>  \param[in]      mean    :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`, representing the center of the distribution.<br>
    !>                              (**optional**, default = `[(0., i = 1, size(rand))]`.)
    !>  \param[in]      chol    :   The input `contiguous` matrix of shape `(ndim, ndim)` whose specified triangular `subset` contains
    !>                              the [Cholesky Factorization](@ref pm_matrixChol) of the Gramian matrix of the MVUE distribution.<br>
    !>                              (**optional**, the default is the Identity matrix of rank `ndim`. It must be present <b>if and only if</b> the input argument `subset` is also present.)
    !>  \param[in]      subset  :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                          implying that the upper-diagonal triangular block of the input `chol` must be used while the lower subset is not referenced.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                          implying that the lower-diagonal triangular block of the input `chol` must be used while the upper subset is not referenced.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `chol` is present.)
    !>
    !>  \interface{setUnifEllRand}
    !>  \code{.F90}
    !>
    !>      use pm_distUnifEll, only: setUnifEllRand
    !>
    !>      ! single vector, using default rng
    !>
    !>      call setUnifEllRand(rand(1:ndim), mean(1:ndim))
    !>      call setUnifEllRand(rand(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifEllRand(rand(1:ndim), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! single vector, using custom rng
    !>
    !>      call setUnifEllRand(rng, rand(1:ndim), mean(1:ndim))
    !>      call setUnifEllRand(rng, rand(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifEllRand(rng, rand(1:ndim), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using default rng
    !>
    !>      call setUnifEllRand(rand(1:ndim, 1:nsam), mean(1:ndim))
    !>      call setUnifEllRand(rand(1:ndim, 1:nsam), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifEllRand(rand(1:ndim, 1:nsam), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using custom rng
    !>
    !>      call setUnifEllRand(rng, rand(1:ndim, 1:nsam), mean(1:ndim))
    !>      call setUnifEllRand(rng, rand(1:ndim, 1:nsam), chol(1:ndim, 1:ndim), subset)
    !>      call setUnifEllRand(rng, rand(1:ndim, 1:nsam), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mean, 1) == size(rand, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(chol) == size(rand, 1))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>  The procedures of this generic interface are `pure` when the input argument `rng` is set to
    !>  [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type) and the compile-time macro `CHECK_ENABLED` is set to `0` or is undefined.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getMultiNormLogPDF](@ref pm_distMultiNorm::getMultiNormLogPDF)<br>
    !>
    !>  \example{setUnifEllRand}
    !>  \include{lineno} example/pm_distUnifEll/setUnifEllRand/main.F90
    !>  \compilef{setUnifEllRand}
    !>  \output{setUnifEllRand}
    !>  \include{lineno} example/pm_distUnifEll/setUnifEllRand/main.out.F90
    !>  \postproc{setUnifEllRand}
    !>  \include{lineno} example/pm_distUnifEll/setUnifEllRand/main.py
    !>  \vis{setUnifEllRand}
    !>  \image html pm_distUnifEll/setUnifEllRand/setUnifEllRand.RK.png width=700
    !>  \image html pm_distUnifEll/setUnifEllRand/setUnifEllRandMean.RK.png width=700
    !>  \image html pm_distUnifEll/setUnifEllRand/setUnifEllRandChol.RK.png width=700
    !>  \image html pm_distUnifEll/setUnifEllRand/setUnifEllRandMeanChol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distUnifEll](@ref test_pm_distUnifEll)
    !>
    !>  \final{setUnifEllRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D1 RNGD

    interface setUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D1_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D1_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D1_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D1_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D1_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D1_RK5(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D1_RK4(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D1_RK3(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D1_RK2(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D1_RK1(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D1_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D1_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D1_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D1_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D1_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D1_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D1_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D1_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D1_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D1_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D1_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D1_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D1_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D1_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D1_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D1_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D1_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D1_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D1_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D1_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGF

    interface setUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D1_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D1_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D1_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D1_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D1_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1 RNGX

    interface setUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D1_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D1_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D1_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D1_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D1_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGD

    interface setUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D2_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D2_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D2_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D2_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_DM_DC_XXX_D2_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D2_RK5(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D2_RK4(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D2_RK3(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D2_RK2(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_AM_DC_XXX_D2_RK1(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D2_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D2_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D2_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D2_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_UXD_D2_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D2_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D2_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D2_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D2_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_DM_AC_XLD_D2_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D2_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D2_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D2_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D2_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_UXD_D2_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D2_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D2_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D2_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D2_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGD_AM_AC_XLD_D2_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGD_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGF

    interface setUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_DM_DC_XXX_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D2_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D2_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D2_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D2_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_AM_DC_XXX_D2_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_UXD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_DM_AC_XLD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_UXD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMUR_RNGF_AM_AC_XLD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGF_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D2 RNGX

    interface setUnifEllRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_DM_DC_XXX_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D2_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D2_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D2_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D2_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_AM_DC_XXX_D2_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_UXD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_DM_AC_XLD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_UXD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMUR_RNGX_AM_AC_XLD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMUR_RNGX_AM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distUnifEll