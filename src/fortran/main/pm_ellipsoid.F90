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
!>
!>  \brief
!>  This module contains classes and procedures for setting up and
!>  computing the properties of the hyper-ellipsoids in arbitrary dimensions.
!>
!>  \details
!>  An ellipsoid in Euclidean geomtery \f$\ell\f$ is defined by its representative Gramian matrix \f$\Sigma_\ell\f$,
!>  containing all points in \f$\mathbb{R}^\ndim\f$ that satisfy,
!>  \f{equation}{
!>      \large
!>      (X - \mu_\ell)^T ~ \Sigma_\ell^{-1} ~ (X - \mu_\ell) \leq 1 ~,
!>  \f}
!>  where \f$\mu_\ell\f$ represents the center of the ellipsoid,
!>  \f$(X - \mu_\ell)^T\f$ is the transpose of the vector \f$(X - \mu_\ell)\f$,
!>  and \f$\Sigma_\ell^{-1}\f$ is the inverse of the matrix \f$\Sigma_\ell\f$.<br>
!>
!>  The volume of this ellipsoid is given by,
!>  \f{equation}{
!>      \large
!>      V(\ell) = V_\ndim \sqrt{\left| \Sigma_\ell \right|} ~,
!>  \f}
!>  where \f$\left|\Sigma_\ell\right|\f$ is the determinant of \f$\Sigma_\ell\f$ and,
!>  \f{equation}{
!>      \large
!>      V_\ndim =
!>      \frac{\pi^{\ndim / 2}}{\up\Gamma(1 + \ndim / 2)} =
!>      \begin{cases}
!>          \frac{1}{(\ndim/2)!} \pi^{\ndim/2} & \text{if $\ndim$ is even} \\\\
!>          2^\ndim \frac{1}{\ndim!} \big( \frac{\ndim-1}{2} \big)! ~ \pi^{(\ndim-1)/2} & \text{if $\ndim$ is odd}
!>      \end{cases}
!>  \f}
!>  is the volume of an \f$\ndim\f$-ball (that is, a unit-radius \f$\ndim\f$-dimensional hyper-ball).<br>
!>  It is readily seen that the corresponding unit-volume ellipsoid \f$\widehat\ell\f$ has the representative Gramian matrix,
!>  \f{equation}{
!>      \large
!>      \Sigma_{\widehat\ell} = V_\ell^{-2/\ndim} ~ \Sigma_\ell ~.
!>  \f}
!>  More generally, to scale an ellipsoid \f$\ell\f$ by some factor \f$\alpha\f$ along each coordinate axis,
!>  it suffices to be used the new scaled ellipsoid \f$\ell^*\f$ with the representative Gramian matrix,
!>  \f{equation}{
!>      \large
!>      \Sigma_{\ell^*} = \alpha^2 ~ \Sigma_\ell ~,
!>  \f}
!>  in which case, the volume of \f$\ell^*\f$ becomes,
!>  \f{equation}{
!>      \large
!>      V_{\ell^*} = \alpha^\ndim ~ V_\ell ~.
!>  \f}
!>
!>  <b>The surface area of an \f$\ndim\f$-ball</b><br>
!>  The surface area \f$S\f$ of an \f$\ndim\f$-dimensional ball of unit-radius is related to its volume \f$V\f$ as,<br>
!>  \f{equation}{
!>      S_\ndim = \ndim \times V_\ndim ~,
!>  \f}
!>  where the natural logarithm of \f$V\f$ is returned by
!>  [getLogVolUnitBall](@ref pm_ellipsoid::getLogVolUnitBall) and
!>  [setLogVolUnitBall](@ref pm_ellipsoid::setLogVolUnitBall).<br>
!>
!>  \note
!>  <ol>
!>      <li>    An ellipsoid with a diagonal representative Gramian matrix \f$\Sigma_\ell\f$ is an
!>              uncorrelated **hyper-ellipsoid** whose axes are parallel to the coordinate axes.<br>
!>      <li>    An uncorrelated hyper-ellipsoid with equal diagonal elements in its
!>              representative Gramian matrix \f$\Sigma_\ell\f$ is a **hyper-ball**.<br>
!>  </ol>
!>
!>  \see
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  [n-sphere](https://en.wikipedia.org/wiki/N-sphere)<br>
!>  [n-ball](https://en.wikipedia.org/wiki/N-ball)<br>
!>
!>  \test
!>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
!>
!>  \todo
!>  \pvhigh
!>  The excluded procedure `getLogVolUnion()` in this module needs cleanup and merging with this module.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_ellipsoid

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter     :: MODULE_NAME = "@pm_ellipsoid"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the volume of an \f$\ndim\f$-dimensional ball of unit-radius.
    !>
    !>  \details
    !>  This generic functional interface is an exact functional-interface replacement for the
    !>  generic subroutine interface [setVolUnitBall](@ref setVolUnitBall).<br>
    !>
    !>  \param[in]  ndim    :   The input scalar of the same type and kind as the output `volUnitBall`,
    !>                          containing the number of dimensions of the unit-radius hyper-ball.<br>
    !>
    !>  \return
    !>  `volUnitBall`       :   The output scalar (or array of the same rank as other input array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing natural logarithm of the volume of the unit-radius hyper-ball.<br>
    !>
    !>  \interface{getVolUnitBall}
    !>  \code{.F90}
    !>
    !>      use pm_ellipsoid, only: getVolUnitBall
    !>
    !>      volUnitBall = getVolUnitBall(ndim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= ndim` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setVolUnitBall](@ref setVolUnitBall)<br>
    !>
    !>  \example{getVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/getVolUnitBall/main.F90
    !>  \compilef{getVolUnitBall}
    !>  \output{getVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/getVolUnitBall/main.out.F90
    !>  \postproc{getVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/getVolUnitBall/main.py
    !>  \vis{getVolUnitBall}
    !>  \image html example/pm_ellipsoid/getVolUnitBall/getVolUnitBall.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \final{getVolUnitBall}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getVolUnitBall

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getVolUnitBallIter_RK5(ndim) result(volUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolUnitBallIter_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: volUnitBall
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getVolUnitBallIter_RK4(ndim) result(volUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolUnitBallIter_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: volUnitBall
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getVolUnitBallIter_RK3(ndim) result(volUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolUnitBallIter_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: volUnitBall
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getVolUnitBallIter_RK2(ndim) result(volUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolUnitBallIter_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: volUnitBall
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getVolUnitBallIter_RK1(ndim) result(volUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolUnitBallIter_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: volUnitBall
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !>  \brief
    !>  Return the volume of an \f$\ndim\f$-dimensional ball of unit-radius.
    !>
    !>  \details
    !>  The computation of the volume of an n-ball requires involves [factorials](@ref pm_mathFactorial)
    !>  which are computed in the procedures of this generic interface iteratively.<br>
    !>  This generic subroutine interface is an exact functional-interface replacement for the
    !>  generic functional interface [getVolUnitBall](@ref getVolUnitBall).<br>
    !>
    !>  \param[out] volUnitBall     :   The output scalar (or array of the same rank as other input array-like arguments) of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing natural logarithm of the volume of the unit-radius hyper-ball.<br>
    !>  \param[in]  ndim            :   The input scalar containing the number of dimensions of the unit-radius hyper-ball.<br>
    !>                                  It can be,
    !>                                  <ol>
    !>                                      <li>    of type `integer` of default kind \IK.
    !>                                  </ol>
    !>
    !>  \interface{setVolUnitBall}
    !>  \code{.F90}
    !>
    !>      use pm_ellipsoid, only: setVolUnitBall
    !>
    !>      call setVolUnitBall(volUnitBall, ndim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= ndim` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogVolUnitBall](@ref getLogVolUnitBall)<br>
    !>  [Volume of an n-ball](https://en.wikipedia.org/wiki/Volume_of_an_n-ball)<br>
    !>  [Particular values of the Gamma function](https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function)<br>
    !>
    !>  \example{setVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/setVolUnitBall/main.F90
    !>  \compilef{setVolUnitBall}
    !>  \output{setVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/setVolUnitBall/main.out.F90
    !>  \postproc{setVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/setVolUnitBall/main.py
    !>  \vis{setVolUnitBall}
    !>  \image html example/pm_ellipsoid/setVolUnitBall/setVolUnitBall.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \final{setVolUnitBall}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface setVolUnitBall

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setVolUnitBallIter_RK5(volUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolUnitBallIter_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   :: volUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setVolUnitBallIter_RK4(volUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolUnitBallIter_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   :: volUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setVolUnitBallIter_RK3(volUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolUnitBallIter_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   :: volUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setVolUnitBallIter_RK2(volUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolUnitBallIter_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   :: volUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setVolUnitBallIter_RK1(volUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolUnitBallIter_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   :: volUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the volume of an \f$\ndim\f$-dimensional ball of unit-radius.
    !>
    !>  \details
    !>  This generic functional interface is an exact functional-interface replacement for the
    !>  generic subroutine interface [setLogVolUnitBall](@ref setLogVolUnitBall).<br>
    !>
    !>  \param[in]  ndim    :   The input scalar of the same type and kind as the output `logVolUnitBall`,
    !>                          containing the number of dimensions of the unit-radius hyper-ball.<br>
    !>
    !>  \return
    !>  `logVolUnitBall`    :   The output scalar (or array of the same rank as other input array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing natural logarithm of the volume of the unit-radius hyper-ball.<br>
    !>
    !>  \interface{getLogVolUnitBall}
    !>  \code{.F90}
    !>
    !>      use pm_ellipsoid, only: getLogVolUnitBall
    !>
    !>      logVolUnitBall = getLogVolUnitBall(ndim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= ndim` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogVolUnitBall](@ref setLogVolUnitBall)<br>
    !>
    !>  \example{getLogVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/getLogVolUnitBall/main.F90
    !>  \compilef{getLogVolUnitBall}
    !>  \output{getLogVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/getLogVolUnitBall/main.out.F90
    !>  \postproc{getLogVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/getLogVolUnitBall/main.py
    !>  \vis{getLogVolUnitBall}
    !>  \image html example/pm_ellipsoid/getLogVolUnitBall/getLogVolUnitBall.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \final{getLogVolUnitBall}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getLogVolUnitBall

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogVolUnitBallGamm_RK5(ndim) result(logVolUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnitBallGamm_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: logVolUnitBall
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogVolUnitBallGamm_RK4(ndim) result(logVolUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnitBallGamm_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: logVolUnitBall
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogVolUnitBallGamm_RK3(ndim) result(logVolUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnitBallGamm_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: logVolUnitBall
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogVolUnitBallGamm_RK2(ndim) result(logVolUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnitBallGamm_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: logVolUnitBall
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogVolUnitBallGamm_RK1(ndim) result(logVolUnitBall)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnitBallGamm_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    :: ndim
        real(RKG)                           :: logVolUnitBall
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the volume of an \f$\ndim\f$-dimensional ball of unit-radius.
    !>
    !>  \details
    !>  The computation of the volume of an n-ball requires involves [factorials](@ref pm_mathFactorial)
    !>  which are computed in the procedures of this generic interface via the Fortran intrinsic `log_gamma()`.<br>
    !>  This generic subroutine interface is an exact functional-interface replacement for the
    !>  generic functional interface [getLogVolUnitBall](@ref getLogVolUnitBall).<br>
    !>
    !>  \param[out] logVolUnitBall  :   The output scalar (or array of the same rank as other input array-like arguments) of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing natural logarithm of the volume of the unit-radius hyper-ball.<br>
    !>  \param[in]  ndim            :   The input scalar containing the number of dimensions of the unit-radius hyper-ball.<br>
    !>                                  It can be,
    !>                                  <ol>
    !>                                      <li>    of type `integer` of default kind \IK,
    !>                                              in which case the output will be computed using an iterative factorial algorithm, or
    !>                                      <li>    of the same type and kind as the output `logVolUnitBall`,
    !>                                              in which case the output will be computed using the Fortran intrinsic `log_gamma()`.
    !>                                  </ol>
    !>
    !>  \interface{setLogVolUnitBall}
    !>  \code{.F90}
    !>
    !>      use pm_ellipsoid, only: setLogVolUnitBall
    !>
    !>      call setLogVolUnitBall(logVolUnitBall, ndim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= ndim` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogVolUnitBall](@ref getLogVolUnitBall)<br>
    !>  [Volume of an n-ball](https://en.wikipedia.org/wiki/Volume_of_an_n-ball)<br>
    !>  [Particular values of the Gamma function](https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function)<br>
    !>
    !>  \example{setLogVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/setLogVolUnitBall/main.F90
    !>  \compilef{setLogVolUnitBall}
    !>  \output{setLogVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/setLogVolUnitBall/main.out.F90
    !>  \postproc{setLogVolUnitBall}
    !>  \include{lineno} example/pm_ellipsoid/setLogVolUnitBall/main.py
    !>  \vis{setLogVolUnitBall}
    !>  \image html example/pm_ellipsoid/setLogVolUnitBall/setLogVolUnitBall.RK.png width=700
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setLogVolUnitBall, The runtime performance of [setLogVolUnitBall](@ref pm_ellipsoid::setLogVolUnitBall) for `integer` vs. `real` input `ndim`.}
    !>  \include{lineno} benchmark/pm_ellipsoid/setLogVolUnitBall/main.F90
    !>  \compilefb{setLogVolUnitBall}
    !>  \postprocb{setLogVolUnitBall}
    !>  \include{lineno} benchmark/pm_ellipsoid/setLogVolUnitBall/main.py
    !>  \visb{setLogVolUnitBall}
    !>  \image html benchmark/pm_ellipsoid/setLogVolUnitBall/benchmark.setLogVolUnitBall.runtime.png width=1000
    !>  \image html benchmark/pm_ellipsoid/setLogVolUnitBall/benchmark.setLogVolUnitBall.runtime.ratio.png width=1000
    !>  \moralb{setLogVolUnitBall}
    !>      -#  The benchmark procedures named `ndimInt2RK1` and `ndimInt2RK2` call the generic interface [setLogVolUnitBall](@ref pm_ellipsoid::setLogVolUnitBall)
    !>          with a `ndim` argument of type `integer` of kind \IK and return the result as `real` of kind \RKS and \RKD respectively.<br>
    !>          Conversely, the procedures named `ndimReal2RK1` and `ndimReal2RK2` take `ndim` as `real` of kind \RKS and \RKD respectively
    !>          and return results of the same type and kind as `ndim`.<br>
    !>          The first class of procedure interfaces (working with a `ndim` of type `integer`) compute the result using an iterative factorial computation approach.<br>
    !>          The second class of procedure interfaces (working with a `ndim` of type `real`) use the Fortran intrinsic `log_gamma()` to compute the results.<br>
    !>          This performance difference tends to be about a factor of 10 times for scalar `optional` arguments and grows
    !>          substantially to larger factors with switching to increasing-size array-like `optional` dummy arguments.<br>
    !>      -#  From the benchmark results it appears that the interface accepting `ndim` of type \RK32 offers the fastest algorithm, followed by the iterative methods.<br>
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \final{setLogVolUnitBall}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface setLogVolUnitBall

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogVolUnitBallGamm_RK5(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallGamm_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   :: logVolUnitBall
        real(RKG)           , intent(in)    :: ndim
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogVolUnitBallGamm_RK4(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallGamm_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   :: logVolUnitBall
        real(RKG)           , intent(in)    :: ndim
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogVolUnitBallGamm_RK3(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallGamm_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   :: logVolUnitBall
        real(RKG)           , intent(in)    :: ndim
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogVolUnitBallGamm_RK2(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallGamm_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   :: logVolUnitBall
        real(RKG)           , intent(in)    :: ndim
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogVolUnitBallGamm_RK1(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallGamm_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   :: logVolUnitBall
        real(RKG)           , intent(in)    :: ndim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogVolUnitBallIter_RK5(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallIter_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   :: logVolUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogVolUnitBallIter_RK4(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallIter_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   :: logVolUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogVolUnitBallIter_RK3(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallIter_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   :: logVolUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogVolUnitBallIter_RK2(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallIter_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   :: logVolUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogVolUnitBallIter_RK1(logVolUnitBall, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogVolUnitBallIter_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   :: logVolUnitBall
        integer(IK)         , intent(in)    :: ndim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the volume of an \f$\ndim\f$-dimensional ellipsoid.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_ellipsoid](@ref pm_ellipsoid) for computational and algorithmic details.<br>
    !>
    !>  \param[in]  gramian :   The input matrix of the same type and kind as the output `logVolEll`,
    !>                          containing the **upper triangle and diagonal** of the representative
    !>                          Gramian matrix of the ellipsoid.<br>
    !>
    !>  \return
    !>  `logVolEll`         :   The output scalar of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing natural logarithm of the volume of the \f$\ndim\f$-dimensional hyper-ellipsoid.<br>
    !>
    !>  \interface{getLogVolEll}
    !>  \code{.F90}
    !>
    !>      use pm_ellipsoid, only: getLogVolEll
    !>
    !>      logVolEll = getLogVolEll(gramian(1:ndim, 1:ndim))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If the Cholesky factorization of the input Gramian fails,
    !>  the procedures of this generic interface will abort the program by calling `error stop`.<br>
    !>
    !>  \warning
    !>  The condition `size(gramian, 1) == size(gramian, 2)` must hold for the corresponding input arguments.<br>
    !>  The condition \f$0. < \left|\ms{gramian}\right|\f$ must hold for the corresponding input arguments.<br>
    !>  In other words, the input Gramian must be a positive definite matrix.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  Computing the volume of an ellipsoid using its Gramian in a fixed dimension as
    !>  implemented by the procedure of this generic interface are computationally costly.<br>
    !>  The unnecessary costs can be eliminated by precomputing the natural logarithm of the volume of the unit
    !>  ball in the desired dimension once via [setLogVolUnitBall](@ref pm_ellipsoid::setLogVolUnitBall)
    !>  and adding to it the sum of the natural logarithms of the diagonal elements of the Cholesky factorization of
    !>  the representative Gramian matrix of the ellipsoid.<br>
    !>
    !>  \see
    !>  [getLogVolEll](@ref getLogVolEll)<br>
    !>
    !>  \example{getLogVolEll}
    !>  \include{lineno} example/pm_ellipsoid/getLogVolEll/main.F90
    !>  \compilef{getLogVolEll}
    !>  \output{getLogVolEll}
    !>  \include{lineno} example/pm_ellipsoid/getLogVolEll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \todo
    !>  \phigh
    !>  A positive-definiteness runtime check for the `gramian` input argument of this generic interface must be added.<br>
    !>
    !>  \final{getLogVolEll}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getLogVolEll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getLogVolEll_RK5(gramian) result(logVolEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEll_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: gramian(:,:)
        real(RKG)                                   :: logVolEll
    end function
#endif

#if RK4_ENABLED
    PURE module function getLogVolEll_RK4(gramian) result(logVolEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEll_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: gramian(:,:)
        real(RKG)                                   :: logVolEll
    end function
#endif

#if RK3_ENABLED
    PURE module function getLogVolEll_RK3(gramian) result(logVolEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEll_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: gramian(:,:)
        real(RKG)                                   :: logVolEll
    end function
#endif

#if RK2_ENABLED
    PURE module function getLogVolEll_RK2(gramian) result(logVolEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEll_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: gramian(:,:)
        real(RKG)                                   :: logVolEll
    end function
#endif

#if RK1_ENABLED
    PURE module function getLogVolEll_RK1(gramian) result(logVolEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEll_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: gramian(:,:)
        real(RKG)                                   :: logVolEll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` <b>if and only if</b> the input point is a
    !>  member (i.e., inside) of the specified \f$\ndim\f$-dimensional ellipsoid.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_ellipsoid](@ref pm_ellipsoid) for computational and algorithmic details.<br>
    !>
    !>  \param[in]  radius  :   The input positive scalar of the same type and kind as the input argument `point`,
    !>                          containing the full inverse representative Gramian matrix of the \f$\ndim\f$-dimensional ellipsoid.<br>
    !>                          (**optional**. If present, the ellipsoid is assumed to be a ball.
    !>                          It must be present **if and only if** the input argument `invGram` is missing.)
    !>  \param[in]  invGram :   The input square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the input argument `point`,
    !>                          containing the full inverse representative Gramian matrix of the \f$\ndim\f$-dimensional ellipsoid.<br>
    !>                          (**optional**. It must be present **if and only if** the input argument `radius` is missing.)
    !>  \param[in]  center  :   The input vector of shape `(1:ndim)` of the same type and kind as the input `invGram`,
    !>                          containing the center of the ellipsoid.<br>
    !>                          (**optional**. default = `[(0., i = 1, ndim)]`.)
    !>  \param[in]  point   :   The input vector of shape `(1:ndim)` of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the coordinates of a point whose presence within the ellipsoid is to be checked.<br>
    !>
    !>  \return
    !>  `isMember`          :   The output scalar of type `logical` of default kind \LK that is `.true.` <b>if and only if</b>
    !>                          the input point is inside of the specified \f$\ndim\f$-dimensional ellipsoid.<br>
    !>
    !>  \interface{isMemberEll}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_ellipsoid, only: isMemberEll
    !>      logical(LK) :: isMember
    !>
    !>      isMember = isMemberEll(radius, point(1:ndim))
    !>      isMember = isMemberEll(radius, center(1:ndim), point(1:ndim))
    !>      isMember = isMemberEll(invGram(1:ndim, 1:ndim), point(1:ndim)))
    !>      isMember = isMemberEll(invGram(1:ndim, 1:ndim), center(1:ndim), point(1:ndim)))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < radius` must hold for the corresponding input arguments.<br>
    !>  The condition `size(point, 1) == size(center)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(point, 1) == shape(gramian))` must hold for the corresponding input arguments.<br>
    !>  The condition \f$0. < \left|\ms{invGram}\right|\f$ must hold for the corresponding input arguments.<br>
    !>  In other words, the input inverse Gramian must be a positive definite matrix.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [isMemberEll](@ref pm_ellipsoid::isMemberEll)<br>
    !>  [getCountMemberEll](@ref pm_ellipsoid::getCountMemberEll)<br>
    !>
    !>  \example{isMemberEll}
    !>  \include{lineno} example/pm_ellipsoid/isMemberEll/main.F90
    !>  \compilef{isMemberEll}
    !>  \output{isMemberEll}
    !>  \include{lineno} example/pm_ellipsoid/isMemberEll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \todo
    !>  \phigh
    !>  A positive-definiteness runtime check for the `invGram` input argument of this generic interface must be added.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The current implementation of this procedure hard-codes the computation of the Mahalanobis distance,
    !>  whose current implementation relies on the Fortran intrinsic routines.<br>
    !>  This hardcoding should be replaced with potential BLAS alternatives.<br>
    !>
    !>  \final{isMemberEll}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface isMemberEll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMemberSphOrg_RK5(radius, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphOrg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK4_ENABLED
    PURE module function isMemberSphOrg_RK4(radius, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphOrg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK3_ENABLED
    PURE module function isMemberSphOrg_RK3(radius, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphOrg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK2_ENABLED
    PURE module function isMemberSphOrg_RK2(radius, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphOrg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK1_ENABLED
    PURE module function isMemberSphOrg_RK1(radius, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphOrg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMemberSphCen_RK5(radius, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphCen_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK4_ENABLED
    PURE module function isMemberSphCen_RK4(radius, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphCen_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK3_ENABLED
    PURE module function isMemberSphCen_RK3(radius, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphCen_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK2_ENABLED
    PURE module function isMemberSphCen_RK2(radius, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphCen_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

#if RK1_ENABLED
    PURE module function isMemberSphCen_RK1(radius, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberSphCen_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        logical(LK)                                 :: isMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMemberEllOrg_RK5(invGram, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllOrg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK4_ENABLED
    PURE module function isMemberEllOrg_RK4(invGram, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllOrg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK3_ENABLED
    PURE module function isMemberEllOrg_RK3(invGram, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllOrg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK2_ENABLED
    PURE module function isMemberEllOrg_RK2(invGram, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllOrg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK1_ENABLED
    PURE module function isMemberEllOrg_RK1(invGram, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllOrg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMemberEllCen_RK5(invGram, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllCen_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK4_ENABLED
    PURE module function isMemberEllCen_RK4(invGram, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllCen_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK3_ENABLED
    PURE module function isMemberEllCen_RK3(invGram, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllCen_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK2_ENABLED
    PURE module function isMemberEllCen_RK2(invGram, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllCen_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

#if RK1_ENABLED
    PURE module function isMemberEllCen_RK1(invGram, center, point) result(isMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMemberEllCen_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        logical(LK)                                 :: isMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of points that are members (i.e., inside) of the specified \f$\ndim\f$-dimensional ellipsoid.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_ellipsoid](@ref pm_ellipsoid) for computational and algorithmic details.<br>
    !>
    !>  \param[in]  radius  :   The input positive scalar of the same type and kind as the input argument `point`,
    !>                          containing the full inverse representative Gramian matrix of the \f$\ndim\f$-dimensional ellipsoid.<br>
    !>                          (**optional**. If present, the ellipsoid is assumed to be a ball.
    !>                          It must be present **if and only if** the input argument `invGram` is missing.)
    !>  \param[in]  invGram :   The input square matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the input argument `point`,
    !>                          containing the full inverse representative Gramian matrix of the \f$\ndim\f$-dimensional ellipsoid.<br>
    !>                          (**optional**. It must be present **if and only if** the input argument `radius` is missing.)
    !>  \param[in]  center  :   The input vector of shape `(1:ndim)` of the same type and kind as the input `invGram`,
    !>                          containing the center of the ellipsoid.<br>
    !>                          (**optional**. default = `[(0., i = 1, ndim)]`.)
    !>  \param[in]  point   :   The input matrix of shape `(1:ndim, 1:npnt)` of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing `npnt` points whose presence within the ellipsoid is to be checked.<br>
    !>
    !>  \return
    !>  `countMemberEll`    :   The output scalar of type `integer` of default kind \IK,
    !>                          containing the number of points that are within the specified ellipsoid.<br>
    !>
    !>  \interface{getCountMemberEll}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_ellipsoid, only: getCountMemberEll
    !>      integer(IK) :: countMemberEll
    !>
    !>      countMemberEll = getCountMemberEll(radius, point(1:ndim, 1:npnt))
    !>      countMemberEll = getCountMemberEll(radius, center(1:ndim), point(1:ndim, 1:npnt))
    !>      countMemberEll = getCountMemberEll(invGram(1:ndim, 1:ndim), point(1:ndim, 1:npnt)))
    !>      countMemberEll = getCountMemberEll(invGram(1:ndim, 1:ndim), center(1:ndim), point(1:ndim, 1:npnt)))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < radius` must hold for the corresponding input arguments.<br>
    !>  The condition `size(point, 1) == size(center)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(point, 1) == shape(invGram))` must hold for the corresponding input arguments.<br>
    !>  The condition \f$0. < \left|\ms{invGram}\right|\f$ must hold for the corresponding input arguments.<br>
    !>  In other words, the input inverse Gramian must be a positive definite matrix.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [isMemberEll](@ref isMemberEll)<br>
    !>  [getCountMemberEll](@ref getCountMemberEll)<br>
    !>
    !>  \example{getCountMemberEll}
    !>  \include{lineno} example/pm_ellipsoid/getCountMemberEll/main.F90
    !>  \compilef{getCountMemberEll}
    !>  \output{getCountMemberEll}
    !>  \include{lineno} example/pm_ellipsoid/getCountMemberEll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_ellipsoid](@ref test_pm_ellipsoid)<br>
    !>
    !>  \final{getCountMemberEll}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getCountMemberEll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountMemberSphOrg_RK5(radius, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphOrg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountMemberSphOrg_RK4(radius, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphOrg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountMemberSphOrg_RK3(radius, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphOrg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountMemberSphOrg_RK2(radius, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphOrg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountMemberSphOrg_RK1(radius, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphOrg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountMemberSphCen_RK5(radius, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphCen_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountMemberSphCen_RK4(radius, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphCen_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountMemberSphCen_RK3(radius, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphCen_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountMemberSphCen_RK2(radius, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphCen_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountMemberSphCen_RK1(radius, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberSphCen_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)                    :: radius
        integer(IK)                                 :: countMemberEll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountMemberEllOrg_RK5(invGram, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllOrg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountMemberEllOrg_RK4(invGram, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllOrg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountMemberEllOrg_RK3(invGram, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllOrg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountMemberEllOrg_RK2(invGram, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllOrg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountMemberEllOrg_RK1(invGram, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllOrg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getCountMemberEllCen_RK5(invGram, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllCen_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK4_ENABLED
    PURE module function getCountMemberEllCen_RK4(invGram, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllCen_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK3_ENABLED
    PURE module function getCountMemberEllCen_RK3(invGram, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllCen_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK2_ENABLED
    PURE module function getCountMemberEllCen_RK2(invGram, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllCen_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

#if RK1_ENABLED
    PURE module function getCountMemberEllCen_RK1(invGram, center, point) result(countMemberEll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountMemberEllCen_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: point(:,:)
        real(RKG)   , intent(in)    , contiguous    :: center(:)
        real(RKG)   , intent(in)    , contiguous    :: invGram(:,:)
        integer(IK)                                 :: countMemberEll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if 0
!    !>  \brief
!    !>  Generate Compute the *effective* total volume of all bounded ellipsoidal regions together while considering the possible overlaps.
!    !>  This computation is done aggressively via Monte Carlo simulation.
!    !>
!    !> \param[in]       nd              :   The number of dimensions.
!    !> \param[in]       nc              :   The number of input clusters.
!    !> \param[in]       reltol          :   The accuracy tolerance of the computed effective sum of volumes (NOT the log-transformed of the sum).
!    !> \param[in]       logVol          :   An array of size `(nc)` representing the logarithms of the volumes of the ellipsoids normalized to the volume of an nd-ball.
!    !> \param[in]       center          :   An array of size `(nd,nc)` representing the centers of the ellipsoids.
!    !> \param[in]       choDia          :   An array of size `(nd,nc)` representing the diagonal elements of the Cholesky lower triangle of the covariance matrices of the ellipsoids.
!    !> \param[in]       logSumVol       :   `` = sum(logVol)``.
!    !> \param[in]       invCov          :   An array of size `(nd,nd,nc)` representing the inverse covariance matrices of the ellipsoids.
!    !> \param[in]       choLow          :   An array of size `(nd,nd,nc)` whose lower-triangle represents the Cholesky lower triangle of the covariance matrices of the ellipsoids.
!    !> \param[inout]    logSumVolEff    :   The logarithm of the sum of `LogVolNormed`. On output, it will be rewritten with the effective log-sum of the (normalized) volumes.
!    !> \param[in]       ScaleFactorSq   :   An array of size `(nc)` representing the maximum of Mahalanobis-distance-squared of the members of each cluster from its center.
!    !>                                      This is needed if the ellipsoid properties have not been scaled, i.e., `ScaleFactorSq > 1` (**optional**, the default is `ScaleFactorSq(1:nc) = 1`).
!    !
!    !  scratch:
!    !  \param[out]      Overlap         :   An array of size `(nc,nc)` representing the diagonal elements of the Cholesky lower triangle of the covariance matrices of the ellipsoids.
!    subroutine getLogVolUnion   ( reltol & ! LCOV_EXCL_LINE
!                                , logVol & ! LCOV_EXCL_LINE
!                                , center & ! LCOV_EXCL_LINE
!                                , choLow & ! LCOV_EXCL_LINE
!                                , invCov & ! LCOV_EXCL_LINE
!                                , logSumVol & ! LCOV_EXCL_LINE
!                                , logSumVolEff & ! LCOV_EXCL_LINE
!                                , ScaleFactorSq & ! LCOV_EXCL_LINE
!                                !, Overlap & ! LCOV_EXCL_LINE
!                                )
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnion
!#endif
!        use pm_distUnifEll, only: setUnifEllRand, lowDia
!        use pm_mathCumSum, only: setCumSum
!        use pm_kind, only: RK
!        implicit none
!
!        integer(IK) , intent(in)                :: reltol
!        real(RK)    , intent(in)                :: logSumVol
!        real(RK)    , intent(in)                :: logVol(nc)
!        real(RK)    , intent(in)                :: center(nd,nc)
!        real(RK)    , intent(in)                :: invCov(nd,nd,nc)
!        real(RK)    , intent(in)                :: choLow(nd,nd,nc)
!        real(RK)    , intent(out)               :: logSumVolEff
!        real(RK)    , intent(in)    , optional  :: ScaleFactorSq(nc)
!        !integer(IK) , intent(out)   :: Overlap(nc,nc)
!
!        integer(IK)                 :: ic,jc
!        integer(IK)                 :: nsim,isim
!        integer(IK)                 :: membershipCount
!        real(RK)                    :: CumSumVolNormed(nc)
!        real(RK)                    :: ScaleFacSq(nc)
!        real(RK)                    :: VolNormed(nc)
!        real(RK)                    :: isimAccepted
!        real(RK)                    :: simsam(nd)
!        real(RK)                    :: unifrnd
!        real(RK)                    :: mahalSq
!
!        if (present(ScaleFactorSq)) then
!            ScaleFacSq(1:nc) = ScaleFactorSq
!        else
!            ScaleFacSq(1:nc) = 1._RK
!        end if
!
!        ! Compute the expected number of simulated points in each bounded region.
!
!        nsim = nint(1._RK / reltol**2)
!        VolNormed = exp(logVol - logSumVol)
!        call setCumSum(CumSumVolNormed, VolNormed)
!        CumSumVolNormed(1:nc) = CumSumVolNormed(1:nc) / CumSumVolNormed(nc)
!
!#if     CHECK_ENABLED
!        if (abs(CumSumVolNormed(nc) - 1._RK) > 1.e-6) then
!            write(*,*) "abs(CumSumVolNormed(nc) - 1._RK) > 1.e-6", CumSumVolNormed(nc)
!            error stop
!        end if
!#endif
!
!        ! Simulate points within each cluster and compute the effective sum of volumes by excluding overlaps.
!
!        isim = 0._RK
!        isimAccepted = 0._RK
!        loopGeneratePoint: do
!
!            call random_number(unifrnd)
!            ic = minloc(CumSumVolNormed, mask = CumSumVolNormed >= unifrnd, dim = 1)
!
!            call setUnifEllRand( rand = simsam & ! LCOV_EXCL_LINE
!                                    , mean = center(1:nd,ic) & ! LCOV_EXCL_LINE
!                                    , chol = choLow(1:nd,1:nd, ic) & ! LCOV_EXCL_LINE
!                                    , subset = lowDia & ! LCOV_EXCL_LINE
!                                    )
!
!
!            membershipCount = 1_IK
!            do jc = 1, nc
!                if (jc/=ic) then
!                    simsam(1:nd) = simsam(1:nd) - center(1:nd,jc)
!                    mahalSq = dot_product( simsam(1:nd) , matmul(invCov(1:nd,1:nd,jc), simsam(1:nd)) )
!                    simsam(1:nd) = simsam(1:nd) + center(1:nd,jc)
!                    if (mahalSq < ScaleFacSq(jc)) membershipCount = membershipCount + 1
!                end if
!            end do
!
!            isim = isim + 1
!            isimAccepted = isimAccepted + 1._RK / membershipCount
!            if (isim < nsim) cycle loopGeneratePoint
!            exit loopGeneratePoint
!
!        end do loopGeneratePoint
!
!        if (isimAccepted < isim) then
!            logSumVolEff = log(isimAccepted / isim) + logSumVol
!        else
!            logSumVolEff = logSumVol
!        end if
!
!    end subroutine getLogVolUnion
!#endif
!
!
!    !>  \brief
!    !>  Return the (bounding) ellipsoid of the input set of `Point`s.
!    !>
!    !>  \param[in]  ndim          :   The number of dimensions.
!    !>  \param[in]  np          :   The number of input points.
!    !>  \param[in]  Point       :   The input array of points of shape `(ndim,np)`.
!    !>  \param[in]  isBounding  :   A logical input value that if `.true.` will yield the bounding ellipsoid of the input `Point` dataset.
!    !>  \param[out] failed      :   A logical output value that is `.true.` only if the Cholesky factorization fails.
!    !>
!    !>  \return
!    !>  `ellipsoid` : An object of type [ellipsoid_type](@ref ellipsoid_type).
!    !>
!    !>  \warning
!    !>  If any error occurs during the ellipsoid construction, the first element of the `choDia`
!    !>  component of the output `ellipsoid` object will be set to a negative number on return.
!    !>  This is to preserve the purity of the function and to avoid the need for error handling.
!    !>
!    !>
!    !>
!    !>  \remark
!    !>  When `isBounding = .true.`, the `scale`, `scaleSq`, `logVolNormed` are computed as scaled.
!    !>
!    !>  \todo
!    !>  \phigh The `isBounding = .false.` scenario still needs work. Perhaps this argument should be changed to `boundingStatus` with three
!    !>  possible values, corresponding to three different degrees of computations: 1) no bounding, 2) computing scalefactor and logVolNormed,
!    !>  3) computing (2) as well as scaling the properties.
!    function ellipsoid_typer(ndim, np, Point, isBounding, failed) result(ellipsoid)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: ellipsoid_typer
!#endif
!        use pm_kind, only: IK, LK, RK ! LCOV_EXCL_LINE
!        use pm_matrixInv, only: getMatInvFromChoLow
!        use pm_matrixChol, only: setChoLow
!        implicit none
!        integer(IK) , intent(in)            :: ndim, np
!        real(RKG)   , intent(in)            :: Point(ndim,np)
!        logical(LK) , intent(in)            :: isBounding
!        logical(LK) , intent(out)           :: failed
!        type(ellipsoid_type)                :: ellipsoid
!
!        real(RKG)                           :: scaleSqInverse
!        real(RKG)                           :: NormedPoint(ndim,np)
!        integer(IK)                         :: id, ip
!
!        ! Compute the mean.
!
!        ellipsoid%center = sum(Point, dim = 2) / np
!
!        do concurrent(ip = 1:np)
!            NormedPoint(1:ndim,ip) = Point(1:ndim,ip) - ellipsoid%center
!        end do
!
!        ! Compute the upper representative matrix of the cluster representative matrices.
!
!        allocate(ellipsoid%choLowCovUpp(1:ndim,1:ndim), ellipsoid%choDia(1:ndim), source = 0._RK)
!        do ip = 1, np
!            do id = 1, ndim
!                ellipsoid%choLowCovUpp(1:id,id) = ellipsoid%choLowCovUpp(1:id,id) + NormedPoint(1:id,ip) * NormedPoint(id,ip)
!            end do
!        end do
!        ellipsoid%choLowCovUpp(1:ndim,1:ndim) = ellipsoid%choLowCovUpp(1:ndim,1:ndim) / real(np - 1, RK)
!
!        ! Compute the Cholesky Factor of the cluster representative matrices.
!
!        call setChoLow(ellipsoid%choLowCovUpp(1:ndim,1:ndim), ellipsoid%choDia(1:ndim), ndim)
!        failed = ellipsoid%choDia(1) < 0._RK
!        if (failed) return
!
!        ! Compute the inverse of the cluster representative matrices.
!
!        ellipsoid%invCov = getMatInvFromChoLow(ndim, ellipsoid%choLowCovUpp, ellipsoid%choDia)
!
!        ! Compute the mahalSq of the points.
!
!        allocate(ellipsoid%mahalSq(np))
!
!!        if (isBounding) then
!
!            ! Compute the scaleFcator of the bounding region.
!
!            ellipsoid%scaleSq = NEGBIG_RK
!            do ip = 1, np
!                ellipsoid%mahalSq(ip) = dot_product( NormedPoint(1:ndim,ip) , matmul(ellipsoid%invCov, NormedPoint(1:ndim,ip)) )
!                if (ellipsoid%scaleSq < ellipsoid%mahalSq(ip)) ellipsoid%scaleSq = ellipsoid%mahalSq(ip)
!            end do
!
!            ellipsoid%scale = sqrt(ellipsoid%scaleSq)
!            scaleSqInverse = 1._RK / ellipsoid%scaleSq
!
!            ellipsoid%choDia(1:ndim) = ellipsoid%choDia * ellipsoid%scale
!            ellipsoid%invCov(1:ndim,1:ndim) = ellipsoid%invCov * scaleSqInverse
!            do concurrent(id = 1:ndim); ellipsoid%choLowCovUpp(id+1:ndim,id) = ellipsoid%choLowCovUpp(id+1:ndim,id) * ellipsoid%scale; end do
!
!!        else
!!
!!            do concurrent(ip = 1:np)
!!                ellipsoid%mahalSq(ip) = dot_product( NormedPoint(1:ndim,ip) , matmul(ellipsoid%invCov, NormedPoint(1:ndim,ip)) )
!!            end do
!!
!!            ellipsoid%scale = 1._RK
!!            ellipsoid%scaleSq = 1._RK
!!            ellipsoid%logVolNormed = sum(log(ellipsoid%choDia(1:ndim)))
!!
!!        end if
!
!        ellipsoid%logVolNormed = sum(log(ellipsoid%choDia(1:ndim)))
!
!    end function ellipsoid_typer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    pure function getLogLikeFitness(count, logVolRatio) result(logLikeFitness)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogLikeFitness
!#endif
!        use pm_kind, only: RKG
!        implicit none
!        integer(IK) , intent(in)    :: count
!        real(RKG)   , intent(in)    :: logVolRatio
!        real(RKG)                   :: logLikeFitness
!        logLikeFitness = count * (1._RK + logVolRatio - exp(logVolRatio))
!    end function getLogLikeFitness
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_ellipsoid ! LCOV_EXCL_LINE