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
!>  This module contains procedures and generic interfaces for computing the Multivariate Normal Shell density function
!>  or mixtures of such densities with varying parameters.
!>
!>  \details
!>  The NormShell density function is frequently used in testing the efficiency of optimization and sampling algorithms.<br>
!>  Given a value \f$x\in\mathbb{R}^n\f$, the NormShell density function with the <b>(location, scale, scale, scale)</b>
!>  parameters \f$(\bu{\mu}, \bu{\Sigma}, \omega, \rho)\f$ is defined as,
!>  \f{equation}{
!>      \large
!>      f\left( \bu{X} ~|~\bu{\mu}, \bu{\Sigma}, \omega, \rho \right) =
!>      \exp \left( -\frac{\left[\left|(\bu{X} - \bu{\mu})^T ~ \bu{\Sigma}^{-1} ~ (\bu{X} - \bu{\mu})\right| - \rho\right]^2}{2\omega^2} \right) ~,
!>  \f}
!>  where \f$(\bu{\mu}, \bu{\Sigma}, \omega, \rho)\f$ are the **shell center**, **shell covariance matrix**,
!>  **shell width**, and **shell radius** of the distribution respectively.<br>
!>
!>  \see
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  [pm_distEggBox](@ref pm_distEggBox)<br>
!>
!>  \test
!>  [test_pm_distNormShell](@ref test_pm_distNormShell)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distNormShell

    use pm_kind, only: IK, LK, RK, SK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distNormShell"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the NormShell density function value(s) at the specified input point `X`,
    !>  for the specified set of parameters of the single or mixture of NormShell distributions.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distNormShell](@ref pm_distNormShell) for details of the NormShell density function.<br>
    !>
    !>  \param[in]  X               :   The input vector of size `(1:ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the `ndim`-dimensional point at which the function must be evaluated.
    !>  \param[in]  center          :   The input argument of the same type and kind as `X`,
    !>                                  representing the center (location parameter: \f$\mu\f$) of the density function(s).<br>
    !>                                  It can be,
    !>                                  <ol>
    !>                                      <li>    a `contiguous` vector of the same rank and size as the input `X`: `center(1:ndim)`.<br>
    !>                                      <li>    a `contiguous` matrix whose first dimension is of the same size as the input `X`,
    !>                                              and its second dimension is the number of unique density functions in the mixture: `center(1:ndim, 1:nmix)`.
    !>                                  </ol>
    !>                                  (**optional**, default = `0.`. It must be present **if and only if** the input argument `invCov` is also present.)
    !>  \param[in]  invCov          :   The input argument of the same type and kind as `X`,
    !>                                  representing the inverse covariance matrix(es) (\f$\Sigma^{-1}\f$) of the Normal distribution(s) of the density function(s).<br>
    !>                                  It can be,
    !>                                  <ol>
    !>                                      <li>    a `contiguous` square matrix of the same rank as the size of the input `X`: `invCov(1:ndim, 1:ndim)`.<br>
    !>                                      <li>    a `contiguous` array of rank `3` whose first and second dimensions are of the same size as the input `X`,
    !>                                              and its third dimension is the number of unique density functions in the mixture: `invCov(1:ndim, 1:ndim, 1:nmix)`.
    !>                                  </ol>
    !>                                  (**optional**, default = the identity matrix. It must be present **if and only if** the input argument `center` is also present.)
    !>  \param[in]  width           :   The input positive-valued argument of the same type and kind as the input `X`,
    !>                                  representing the width of the shell(s) as specified by the shape parameter (\f$\omega\f$) of the density function(s).<br>
    !>                                  It can be,
    !>                                  <ol>
    !>                                      <li>    a scalar representing the width of the single density function.
    !>                                      <li>    a `contiguous` vector of size `nmix`,
    !>                                              representing the widths of individual density function(s) in the mixture.
    !>                                  </ol>
    !>                                  (**optional**, default = `1.`)
    !>  \param[in]  radius          :   The input positive-valued argument of the same type and kind as the input `X`,
    !>                                  representing the radius of the shell(s) as specified by the shape parameter (\f$\rho\f$) of the density function(s).<br>
    !>                                  It can be,
    !>                                  <ol>
    !>                                      <li>    a scalar representing the radius of the single density function.
    !>                                      <li>    a `contiguous` vector of size `nmix`,
    !>                                              representing the radii of individual density function(s) in the mixture.
    !>                                  </ol>
    !>                                  (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logUDF`                    :   The output of the same type and kind as the input argument `X` representing the natural logarithm
    !>                                  of the value(s) of the density function(s) at the specified location `X` and parameters sets.<br>
    !>                                  It is,
    !>                                  <ol>
    !>                                      <li>    a scalar **if and only if** the dimensionality of the
    !>                                              input arguments corresponds to a single density function.
    !>                                      <li>    a vector of size `nmix` **if and only if** the dimensionality of the
    !>                                              input arguments corresponds to a mixture of `nmix` density functions.
    !>                                  </ol>
    !>
    !>  \interface{getNormShellLogUDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNormShell, only: getNormShellLogUDF
    !>
    !>      ! single density function.
    !>
    !>      logUDF = getNormShellLogUDF(X(1:ndim), width = width, radius = radius)
    !>      logUDF = getNormShellLogUDF(X(1:ndim), center(1:ndim), invCov(1:ndim, 1:ndim), width = width, radius = radius)
    !>
    !>      ! mixture of density functions.
    !>
    !>      logUDF(1:nmix) = getNormShellLogUDF(X(1:ndim), width = width(1:nmix), radius = radius(1:nmix))
    !>      logUDF(1:nmix) = getNormShellLogUDF(X(1:ndim), center(1:ndim, 1:nmix), invCov(1:ndim, 1:ndim, 1:nmix), width = width(1:nmix), radius = radius(1:nmix))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all([0. < width])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([0. < radius])` must hold for the corresponding input arguments.<br>
    !>  The condition `size(X, 1) == size(center, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(X, 1) == [size(invCov, 1), size(invCov, 2)])` must hold for the corresponding input arguments.<br>
    !>  The condition `size(invCov, rank(invCov)) == size(center, rank(invCov))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \example{getNormShellLogUDF}
    !>  \include{lineno} example/pm_distNormShell/getNormShellLogUDF/main.F90
    !>  \compilef{getNormShellLogUDF}
    !>  \output{getNormShellLogUDF}
    !>  \include{lineno} example/pm_distNormShell/getNormShellLogUDF/main.out.F90
    !>  \postproc{getNormShellLogUDF}
    !>  \include{lineno} example/pm_distNormShell/getNormShellLogUDF/main.py
    !>  \vis{getNormShellLogUDF}
    !>  \image html pm_distNormShell/getNormShellLogUDF/getNormShellLogUDF.D1.RK.png width=700
    !>  \image html pm_distNormShell/getNormShellLogUDF/getNormShellLogUDF.D2.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNormShell](@ref test_pm_distNormShell)
    !>
    !>  \final{getNormShellLogUDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNormShellLogUDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormShellLogUDFODD_D1_RK5(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFODD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormShellLogUDFODD_D1_RK4(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFODD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormShellLogUDFODD_D1_RK3(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFODD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormShellLogUDFODD_D1_RK2(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFODD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormShellLogUDFODD_D1_RK1(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFODD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormShellLogUDFOCI_D1_RK5(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFOCI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:), invCov(:,:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormShellLogUDFOCI_D1_RK4(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFOCI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:), invCov(:,:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormShellLogUDFOCI_D1_RK3(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFOCI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:), invCov(:,:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormShellLogUDFOCI_D1_RK2(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFOCI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:), invCov(:,:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormShellLogUDFOCI_D1_RK1(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFOCI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:), invCov(:,:)
        real(RKG)   , intent(in)                , optional  :: width, radius
        real(RKG)                                           :: logUDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormShellLogUDFMDD_D1_RK5(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(width, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormShellLogUDFMDD_D1_RK4(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(width, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormShellLogUDFMDD_D1_RK3(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(width, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormShellLogUDFMDD_D1_RK2(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(width, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormShellLogUDFMDD_D1_RK1(X, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(width, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getNormShellLogUDFMCI_D1_RK5(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMCI_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:,:), invCov(:,:,:)
        real(RKG)   , intent(in)    , contiguous, optional  :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(center, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getNormShellLogUDFMCI_D1_RK4(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMCI_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:,:), invCov(:,:,:)
        real(RKG)   , intent(in)    , contiguous, optional  :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(center, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getNormShellLogUDFMCI_D1_RK3(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMCI_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:,:), invCov(:,:,:)
        real(RKG)   , intent(in)    , contiguous, optional  :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(center, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getNormShellLogUDFMCI_D1_RK2(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMCI_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:,:), invCov(:,:,:)
        real(RKG)   , intent(in)    , contiguous, optional  :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(center, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getNormShellLogUDFMCI_D1_RK1(X, center, invCov, width, radius) result(logUDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormShellLogUDFMCI_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:), center(:,:), invCov(:,:,:)
        real(RKG)   , intent(in)    , contiguous, optional  :: width(:), radius(:)
        real(RKG)                                           :: logUDF(size(center, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distNormShell ! LCOV_EXCL_LINE