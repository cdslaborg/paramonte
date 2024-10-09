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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>MultiVariate Normal (MVN) distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>MultiVariate Normal distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The PDF of the MVN distribution with the mean vector \f$\bu{\mu}\f$ and covariance matrix \f$\bu{\Sigma}\f$
!>  at a given input point \f$X\f$ is defined as,
!>  \f{equation}{
!>      \large
!>      \pi\big(\bu{X} ~|~\bu{\mu}, \bu{\Sigma}\big) =
!>      \frac{1}{\sqrt{\big| 2\pi\bu{\Sigma} \big|}} ~
!>      \exp\bigg( -\frac{1}{2}(\bu{X}-\bu{\mu})^T ~ \bu{\Sigma}^{-1} ~ (\bu{X}-\bu{\mu}) \bigg) ~,
!>  \f}
!>  which is defined <b>if and only if</b> the \f$\bu{\Sigma}\f$ is positive-definite.
!>  The term \f$\large \frac{1}{\sqrt{\big| 2\pi\bu{\Sigma} \big|}}\f$ is
!>  the **Normalization Factor** or the **Normalization Constant** of the MVN PDF whose natural logarithm
!>  can be computed via [getMultiNormLogPDFNF](@ref pm_distMultiNorm::getMultiNormLogPDFNF) in this module.
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [pm_distMultiNorm](@ref pm_distMultiNorm)<br>
!>  [pm_distUnifEll](@ref pm_distUnifEll)<br>
!>  [pm_distUnifPar](@ref pm_distUnifPar)<br>
!>  Applied Multivariate Statistical Analysis, Johnson, Wichern, 1998, 4th ed.<br>
!>
!>  \test
!>  [test_pm_distMultiNorm](@ref test_pm_distMultiNorm)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 20, 2009, 9:12 PM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distMultiNorm

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf_type
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distMultiNorm"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type MultiVariate Normal (MVN)
    !>  as defined in the description of [pm_distMultiNorm](@ref pm_distMultiNorm).
    !>
    !>  \details
    !>  See the documentation of [pm_distMultiNorm](@ref pm_distMultiNorm) for the definition of the MultiVariate Normal (MVN) distribution.
    !>
    !>  \interface{distMultiNorm_type}
    !>  \code{.F90}
    !>
    !>      use pm_distMultiNorm, only: distMultiNorm_type
    !>      type(distMultiNorm_type) :: distMultiNorm
    !>
    !>      distMultiNorm = distMultiNorm_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distMultiNorm](@ref test_pm_distMultiNorm)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distMultiNorm_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distMultiNorm_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the normalization coefficient of the Probability Density Function (PDF)
    !>  of the MultiVariate Normal distribution as defined in the description of [pm_distMultiNorm](@ref pm_distMultiNorm).
    !>
    !>  \details
    !>  See the documentation of [pm_distMultiNorm](@ref pm_distMultiNorm) for the definition of the Normalization Factor of the MVN PDF.
    !>
    !>  \param[in]  ndim                :   The input positive-valued scalar or array of the same shape as other array-like arguments,
    !>                                      of type `integer` of default kind \IK, representing the number of dimensions of the MVN distribution.<br>
    !>                                      (**optional**. It must be present <b>if and only if</b> `invCov` is missing.)
    !>  \param[in]  logSqrtDetInvCov    :   The input scalar or array of the same shape as other array-like arguments,
    !>                                      of the same type and kind as the output `logPDFNF` representing the natural logarithm of the square
    !>                                      root of the determinant of the inverse of the covariance matrix of the target MVN distribution.<br>
    !>                                      This input argument can be obtained by passing the covariance matrix of the MVN distribution to
    !>                                      [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog) or [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog) and negating the result.
    !>                                      \code{.F90}
    !>                                          logSqrtDetInvCov = -getMatDetSqrtLog(cov)
    !>                                      \endcode
    !>                                      (**optional**, default = `0`)
    !>  \param[in]  invCov              :   The input square matrix (of shape `(ndim,ndim)`) of the same type and kind as the output `logPDFNF`
    !>                                      containing the inverse of the covariance matrix of the target MVN distribution.<br>
    !>                                      (**optional**. It must be present <b>if and only if</b> `ndim` is missing.)
    !>  \param[out] info                :   The output scalar `integer` of default kind \IK. On output, it is `0` <b>if and only if</b> the computation
    !>                                      of the determinant of the inverse of the covariance matrix of the target MVN distribution succeeds.<br>
    !>                                      Otherwise, it is set to the order of the leading minor of the specified input subset of `mat` that is not positive definite,
    !>                                      indicating the occurrence of an error.<br>
    !>                                      (**optional**. It can be present only if `invCov` input argument is also present.
    !>                                      If missing and the Cholesky factorization fails, the program will halt by calling `error stop`.)
    !>
    !>  \return
    !>  `logPDFNF`                      :   The output scalar or array of the same shape as the input array-like arguments, of
    !>                                      <ol>
    !>                                          <li>    type `real` of kind \RKALL,
    !>                                      </ol>
    !>                                      containing the normalization factor of the MVN distribution.
    !>
    !>  \interface{getMultiNormLogPDFNF}
    !>  \code{.F90}
    !>
    !>      use pm_distMultiNorm, only: getMultiNormLogPDFNF
    !>
    !>      logPDFNF = getMultiNormLogPDFNF(invCov)
    !>      logPDFNF = getMultiNormLogPDFNF(invCov, info)
    !>      logPDFNF = getMultiNormLogPDFNF(ndim, logSqrtDetInvCov)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ndim` must hold.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the output optional argument `info` is present.<br>
    !>
    !>  \elemental
    !>  The procedures of this generic interface are always non-`elemental` when the input optional argument `invCov` is present.<br>
    !>
    !>  \example{getMultiNormLogPDFNF}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormLogPDFNF/main.F90
    !>  \compilef{getMultiNormLogPDFNF}
    !>  \output{getMultiNormLogPDFNF}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormLogPDFNF/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distMultiNorm](@ref test_pm_distMultiNorm)
    !>
    !>  \final{getMultiNormLogPDFNF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getMultiNormLogPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFNFI_RK5(invCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFNFI_RK4(invCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFNFI_RK3(invCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFNFI_RK2(invCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFNFI_RK1(invCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        real(RKG)                               :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMultiNormLogPDFNFIF_RK5(invCov, info) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFIF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        integer(IK) , intent(out)               :: info
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    impure module function getMultiNormLogPDFNFIF_RK4(invCov, info) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFIF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        integer(IK) , intent(out)               :: info
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    impure module function getMultiNormLogPDFNFIF_RK3(invCov, info) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFIF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        integer(IK) , intent(out)               :: info
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    impure module function getMultiNormLogPDFNFIF_RK2(invCov, info) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFIF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        integer(IK) , intent(out)               :: info
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    impure module function getMultiNormLogPDFNFIF_RK1(invCov, info) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFIF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: invCov(:,:)
        integer(IK) , intent(out)               :: info
        real(RKG)                               :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getMultiNormLogPDFNFND_RK5(ndim, logSqrtDetInvCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFND_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                :: ndim
        real(RKG)   , intent(in)                :: logSqrtDetInvCov
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getMultiNormLogPDFNFND_RK4(ndim, logSqrtDetInvCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFND_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                :: ndim
        real(RKG)   , intent(in)                :: logSqrtDetInvCov
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getMultiNormLogPDFNFND_RK3(ndim, logSqrtDetInvCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFND_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                :: ndim
        real(RKG)   , intent(in)                :: logSqrtDetInvCov
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getMultiNormLogPDFNFND_RK2(ndim, logSqrtDetInvCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFND_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                :: ndim
        real(RKG)   , intent(in)                :: logSqrtDetInvCov
        real(RKG)                               :: logPDFNF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getMultiNormLogPDFNFND_RK1(ndim, logSqrtDetInvCov) result(logPDFNF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFNFND_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                :: ndim
        real(RKG)   , intent(in)                :: logSqrtDetInvCov
        real(RKG)                               :: logPDFNF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF)
    !>  of the MultiVariate Normal distribution as defined in the description of [pm_distMultiNorm](@ref pm_distMultiNorm).
    !>
    !>  \details
    !>  See the documentation of [pm_distMultiNorm](@ref pm_distMultiNorm) for the definition of the MVN PDF.
    !>
    !>  \param[in]  X               :   The input array of shape `(:)` of type `real` of kind \RKALL of arbitrary size
    !>                                  (`ndim` : representing the number of dimensions of the MVN) at which the log(PDF)
    !>                                  of the MVN must be computed.
    !>  \param[in]  mean            :   The input array of the same type kind, and shape as the input `X` containing the mean of the MVN distribution.<br>
    !>                                  (**optional**, default = `[( 0., integer :: i = 1, size(X) )]`)
    !>  \param[in]  invCov          :   The input square matrix (of shape `(ndim,ndim)`) of the same type and kind as the input `X`,
    !>                                  containing the inverse of the covariance matrix of the MVN distribution.<br>
    !>                                  (**optional**, default = Identity Matrix : [getMatInit([size(X), size(X)], uppLowDia_type(), 0., 0., 1.)](@ref pm_matrixInit::getMatInit))
    !>  \param[in]  logPDFNF        :   The input scalar of the same type and kind as the input `X`
    !>                                  containing the normalization factor of the MVN distribution.<br>
    !>                                  Specifying this input argument can lead to faster runtime when the log(PDF) must be computed
    !>                                  for multiple different `X` values given the same inverse covariance matrix for the MVN PDF.<br>
    !>                                  It can be obtained by calling [getMultiNormLogPDFNF](@ref pm_distMultiNorm::getMultiNormLogPDFNF).<br>
    !>                                  (**optional**, default = `getMultiNormLogPDFNF(ndim = size(X), logSqrtDetInvCov)`)
    !>
    !>  \return
    !>  `logPDF`                    :   The output scalar or array of the same shape as the input array-like arguments, of the same type
    !>                                  and kind as the input `logDetInvCovMat` representing the normalization factor of the MVN distribution.
    !>
    !>  \interface{getMultiNormLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distMultiNorm, only: getMultiNormLogPDF
    !>
    !>      logPDF = getMultiNormLogPDF(X(1:ndim))
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), mean(1:ndim))
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), invCov(1:ndim, 1:ndim))
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), mean(1:ndim), invCov(1:ndim, 1:ndim))
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), mean(1:ndim), invCov(1:ndim, 1:ndim), logPDFNF)
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), invCov(1:ndim, 1:ndim), logPDFNF)
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), mean(1:ndim), logPDFNF)
    !>      logPDF = getMultiNormLogPDF(X(1:ndim), logPDFNF)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(size(mean) == size(X, 1))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(invCov) == size(X, 1))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getMultiNormLogPDFNF](@ref pm_distMultiNorm::getMultiNormLogPDFNF)<br>
    !>
    !>
    !>  \example{getMultiNormLogPDF}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormLogPDF/main.F90
    !>  \compilef{getMultiNormLogPDF}
    !>  \output{getMultiNormLogPDF}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormLogPDF/main.out.F90
    !>  \postproc{getMultiNormLogPDF}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormLogPDF/main.py
    !>  \vis{getMultiNormLogPDF}
    !>  \image html pm_distMultiNorm/getMultiNormLogPDF/getMultiNormLogPDF.D2.RK.png width=700
    !>  \test
    !>  [test_pm_distMultiNorm](@ref test_pm_distMultiNorm)<br>
    !>
    !>  \todo
    !>  \pmed
    !>  When the input argument `logPDFNF` is missing and `invCov` is present, there is a
    !>  possibility that the Cholesky Factorization in the computation of `logSqrtDetInvCov` fails.<br>
    !>  In such cases, the procedure will halt the program by calling `error stop`.<br>
    !>  An optional `info` output argument must be added in the future to handle such runtime failures gracefully.
    !>
    !>  \final{getMultiNormLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getMultiNormLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFDDD_D1_RK5(X) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFDDD_D1_RK4(X) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFDDD_D1_RK3(X) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFDDD_D1_RK2(X) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFDDD_D1_RK1(X) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFDDN_D1_RK5(X, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDN_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFDDN_D1_RK4(X, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDN_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFDDN_D1_RK3(X, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDN_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFDDN_D1_RK2(X, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDN_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFDDN_D1_RK1(X, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDDN_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFMDD_D1_RK5(X, mean) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFMDD_D1_RK4(X, mean) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFMDD_D1_RK3(X, mean) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFMDD_D1_RK2(X, mean) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFMDD_D1_RK1(X, mean) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFMDN_D1_RK5(X, mean, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDN_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFMDN_D1_RK4(X, mean, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDN_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFMDN_D1_RK3(X, mean, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDN_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFMDN_D1_RK2(X, mean, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDN_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFMDN_D1_RK1(X, mean, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMDN_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFDID_D1_RK5(X, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFDID_D1_RK4(X, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFDID_D1_RK3(X, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFDID_D1_RK2(X, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFDID_D1_RK1(X, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFDIN_D1_RK5(X, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDIN_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFDIN_D1_RK4(X, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDIN_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFDIN_D1_RK3(X, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDIN_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFDIN_D1_RK2(X, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDIN_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFDIN_D1_RK1(X, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFDIN_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFMID_D1_RK5(X, mean, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMID_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFMID_D1_RK4(X, mean, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMID_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFMID_D1_RK3(X, mean, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMID_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFMID_D1_RK2(X, mean, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMID_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFMID_D1_RK1(X, mean, invCov) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMID_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMultiNormLogPDFMIN_D1_RK5(X, mean, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMIN_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE module function getMultiNormLogPDFMIN_D1_RK4(X, mean, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMIN_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE module function getMultiNormLogPDFMIN_D1_RK3(X, mean, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMIN_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE module function getMultiNormLogPDFMIN_D1_RK2(X, mean, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMIN_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE module function getMultiNormLogPDFMIN_D1_RK1(X, mean, invCov, logPDFNF) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMultiNormLogPDFMIN_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous            :: X(:)
        real(RKG)   , intent(in)    , contiguous            :: mean(:)
        real(RKG)   , intent(in)    , contiguous            :: invCov(:,:)
        real(RKG)   , intent(in)                            :: logPDFNF
        real(RKG)                                           :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a (collection) of random vector(s) of size `ndim` from the `ndim`-dimensional MultiVariate Normal (**MVN**) distribution,
    !>  optionally with the specified input `mean(1:ndim)` and the specified `subset` of the Cholesky Factorization of the Covariance matrix of the MVN distribution.
    !>
    !>  \details
    !>  The procedures of this generic interface are merely wrappers around the subroutine interface [setMultiNormRand](@ref pm_distMultiNorm::setMultiNormRand).<br>
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>                              (**optional**, default = [rngf_type](@ref pm_distUnif::rngf_type).)
    !>  \param[in]      mean    :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`, representing the mean of the Multivariate Normal distribution.<br>
    !>                              (**optional**, default = `[(0., i = 1, size(rand))]`. It must be present if the input argument `chol` is missing.)
    !>  \param[in]      chol    :   The input `contiguous` matrix of shape `(ndim, ndim)` whose specified triangular `subset` contains the [Cholesky Factorization](@ref pm_matrixChol) of the covariance matrix of the MVN distribution.<br>
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
    !>  \param[in]      nsam    :   The input scalar `integer` of default kind \IK containing the number of random MVN vectors to generate.<br>
    !>                              (**optional**. If present, the output `rand` is of rank `2`, otherwise is of rank `1`.)
    !>
    !>  \return
    !>  `rand`                  :   The output vector of
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the MVN distributed random output vector(s):<br>
    !>                              <ol>
    !>                                  <li>    If the input argument `nsam` is missing, then `rand` shall be of shape `(1:ndim)`.<br>
    !>                                  <li>    If the input argument `nsam` is present, then `rand` shall be of shape `(1:ndim, 1:nsam)`.<br>
    !>                              </ol>
    !>
    !>  \interface{getMultiNormRand}
    !>  \code{.F90}
    !>
    !>      use pm_distMultiNorm, only: getMultiNormRand
    !>
    !>      ! single vector, using default rng
    !>
    !>      rand(1:ndim) = getMultiNormRand(mean(1:ndim))
    !>      rand(1:ndim) = getMultiNormRand(chol(1:ndim, 1:ndim), subset)
    !>      rand(1:ndim) = getMultiNormRand(mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! single vector, using custom rng
    !>
    !>      rand(1:ndim) = getMultiNormRand(rng, mean(1:ndim))
    !>      rand(1:ndim) = getMultiNormRand(rng, chol(1:ndim, 1:ndim), subset)
    !>      rand(1:ndim) = getMultiNormRand(rng, mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using default rng
    !>
    !>      rand(1:ndim, 1:nsam) = getMultiNormRand(mean(1:ndim), nsam)
    !>      rand(1:ndim, 1:nsam) = getMultiNormRand(chol(1:ndim, 1:ndim), subset, nsam)
    !>      rand(1:ndim, 1:nsam) = getMultiNormRand(mean(1:ndim), chol(1:ndim, 1:ndim), subset, nsam)
    !>
    !>      ! collection of `nsam` vectors, using custom rng
    !>
    !>      rand(1:ndim, 1:nsam) = getMultiNormRand(rng, mean(1:ndim), nsam)
    !>      rand(1:ndim, 1:nsam) = getMultiNormRand(rng, chol(1:ndim, 1:ndim), subset, nsam)
    !>      rand(1:ndim, 1:nsam) = getMultiNormRand(rng, mean(1:ndim), chol(1:ndim, 1:ndim), subset, nsam)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>
    !>  \example{getMultiNormRand}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormRand/main.F90
    !>  \compilef{getMultiNormRand}
    !>  \output{getMultiNormRand}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormRand/main.out.F90
    !>  \postproc{getMultiNormRand}
    !>  \include{lineno} example/pm_distMultiNorm/getMultiNormRand/main.py
    !>  \vis{getMultiNormRand}
    !>  \image html pm_distMultiNorm/getMultiNormRand/getMultiNormRandMean.RK.png width=700
    !>  \image html pm_distMultiNorm/getMultiNormRand/getMultiNormRandChol.RK.png width=700
    !>  \image html pm_distMultiNorm/getMultiNormRand/getMultiNormRandMeanChol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distMultiNorm](@ref test_pm_distMultiNorm)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      getMNR_RNGD_DM_DC_XXX_D1_RK5()
    !>             |||| || || ||| || |||
    !>             |||| || || ||| || |||
    !>             |||| || || ||| || |||
    !>             |||| || || ||| || |The Kind of the output array.
    !>             |||| || || ||| || The type of the output array: R => Real.
    !>             |||| || || ||| The Dimension of the output array.
    !>             |||| || || The subset of the Cholesky factor: D/U/L => Default/Upper/Lower.
    !>             |||| || The Cholesky factor of the covariance matrix of the distribution: DC/AC => Default/Arbitrary Cholesky
    !>             |||| The mean of the distribution: DM/AM => Default/Arbitrary Mean
    !>             The random number generator: RNG D/F/X => Default/Fortran/Xoroshiro256++
    !>  \endcode
    !>
    !>  \final{getMultiNormRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D1 RNGD

    interface getMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D1_RK5(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D1_RK4(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D1_RK3(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D1_RK2(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D1_RK1(mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D1_RK5(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D1_RK4(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D1_RK3(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D1_RK2(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D1_RK1(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D1_RK5(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D1_RK4(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D1_RK3(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D1_RK2(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D1_RK1(chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D1_RK5(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D1_RK4(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D1_RK3(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D1_RK2(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D1_RK1(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D1_RK5(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D1_RK4(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D1_RK3(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D1_RK2(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D1_RK1(mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D1_RK1
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

    interface getMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D1_RK5(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D1_RK4(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D1_RK3(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D1_RK2(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D1_RK1(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGF_DM_AC_UXD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGF_DM_AC_UXD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGF_DM_AC_UXD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGF_DM_AC_UXD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGF_DM_AC_UXD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D1_RK1
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
    impure module function getMNR_RNGF_DM_AC_XLD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGF_DM_AC_XLD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGF_DM_AC_XLD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGF_DM_AC_XLD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGF_DM_AC_XLD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D1_RK1
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
    impure module function getMNR_RNGF_AM_AC_UXD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGF_AM_AC_UXD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGF_AM_AC_UXD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGF_AM_AC_UXD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGF_AM_AC_UXD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D1_RK1
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
    impure module function getMNR_RNGF_AM_AC_XLD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGF_AM_AC_XLD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGF_AM_AC_XLD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGF_AM_AC_XLD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGF_AM_AC_XLD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D1_RK1
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

    interface getMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D1_RK5(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D1_RK4(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D1_RK3(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D1_RK2(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D1_RK1(rng, mean) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGX_DM_AC_UXD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGX_DM_AC_UXD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGX_DM_AC_UXD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGX_DM_AC_UXD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGX_DM_AC_UXD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D1_RK1
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
    impure module function getMNR_RNGX_DM_AC_XLD_D1_RK5(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGX_DM_AC_XLD_D1_RK4(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGX_DM_AC_XLD_D1_RK3(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGX_DM_AC_XLD_D1_RK2(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGX_DM_AC_XLD_D1_RK1(rng, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D1_RK1
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
    impure module function getMNR_RNGX_AM_AC_UXD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGX_AM_AC_UXD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGX_AM_AC_UXD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGX_AM_AC_UXD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGX_AM_AC_UXD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D1_RK1
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
    impure module function getMNR_RNGX_AM_AC_XLD_D1_RK5(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGX_AM_AC_XLD_D1_RK4(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGX_AM_AC_XLD_D1_RK3(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGX_AM_AC_XLD_D1_RK2(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGX_AM_AC_XLD_D1_RK1(rng, mean, chol, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D1_RK1
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

    interface getMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D2_RK5(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D2_RK4(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D2_RK3(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D2_RK2(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_AM_DC_XXX_D2_RK1(mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D2_RK5(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D2_RK4(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D2_RK3(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D2_RK2(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_DM_AC_UXD_D2_RK1(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_UXD_D2_RK1
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
    impure module function getMNR_RNGD_DM_AC_XLD_D2_RK5(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D2_RK4(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D2_RK3(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D2_RK2(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_DM_AC_XLD_D2_RK1(chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_DM_AC_XLD_D2_RK1
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
    impure module function getMNR_RNGD_AM_AC_UXD_D2_RK5(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D2_RK4(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D2_RK3(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D2_RK2(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_AM_AC_UXD_D2_RK1(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_UXD_D2_RK1
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
    impure module function getMNR_RNGD_AM_AC_XLD_D2_RK5(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D2_RK4(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D2_RK3(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D2_RK2(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)             , intent(in)                    :: nsam
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        real(RKG)                                               :: rand(size(chol, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGD_AM_AC_XLD_D2_RK1(mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGD_AM_AC_XLD_D2_RK1
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

    interface getMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D2_RK5(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D2_RK4(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D2_RK3(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D2_RK2(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGF_AM_DC_XXX_D2_RK1(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_DC_XXX_D2_RK1
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
    impure module function getMNR_RNGF_DM_AC_UXD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D2_RK5
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
    impure module function getMNR_RNGF_DM_AC_UXD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D2_RK4
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
    impure module function getMNR_RNGF_DM_AC_UXD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D2_RK3
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
    impure module function getMNR_RNGF_DM_AC_UXD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D2_RK2
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
    impure module function getMNR_RNGF_DM_AC_UXD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_UXD_D2_RK1
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
    impure module function getMNR_RNGF_DM_AC_XLD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D2_RK5
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
    impure module function getMNR_RNGF_DM_AC_XLD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D2_RK4
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
    impure module function getMNR_RNGF_DM_AC_XLD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D2_RK3
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
    impure module function getMNR_RNGF_DM_AC_XLD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D2_RK2
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
    impure module function getMNR_RNGF_DM_AC_XLD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_DM_AC_XLD_D2_RK1
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
    impure module function getMNR_RNGF_AM_AC_UXD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D2_RK5
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
    impure module function getMNR_RNGF_AM_AC_UXD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D2_RK4
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
    impure module function getMNR_RNGF_AM_AC_UXD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D2_RK3
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
    impure module function getMNR_RNGF_AM_AC_UXD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D2_RK2
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
    impure module function getMNR_RNGF_AM_AC_UXD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_UXD_D2_RK1
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
    impure module function getMNR_RNGF_AM_AC_XLD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D2_RK5
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
    impure module function getMNR_RNGF_AM_AC_XLD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D2_RK4
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
    impure module function getMNR_RNGF_AM_AC_XLD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D2_RK3
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
    impure module function getMNR_RNGF_AM_AC_XLD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D2_RK2
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
    impure module function getMNR_RNGF_AM_AC_XLD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGF_AM_AC_XLD_D2_RK1
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

    interface getMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D2_RK5(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK4_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D2_RK4(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK3_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D2_RK3(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK2_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D2_RK2(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: nsam
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
        real(RKG)                                               :: rand(size(mean, 1, IK), nsam)
    end function
#endif

#if RK1_ENABLED
    impure module function getMNR_RNGX_AM_DC_XXX_D2_RK1(rng, mean, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_DC_XXX_D2_RK1
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
    impure module function getMNR_RNGX_DM_AC_UXD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D2_RK5
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
    impure module function getMNR_RNGX_DM_AC_UXD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D2_RK4
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
    impure module function getMNR_RNGX_DM_AC_UXD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D2_RK3
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
    impure module function getMNR_RNGX_DM_AC_UXD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D2_RK2
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
    impure module function getMNR_RNGX_DM_AC_UXD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_UXD_D2_RK1
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
    impure module function getMNR_RNGX_DM_AC_XLD_D2_RK5(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D2_RK5
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
    impure module function getMNR_RNGX_DM_AC_XLD_D2_RK4(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D2_RK4
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
    impure module function getMNR_RNGX_DM_AC_XLD_D2_RK3(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D2_RK3
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
    impure module function getMNR_RNGX_DM_AC_XLD_D2_RK2(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D2_RK2
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
    impure module function getMNR_RNGX_DM_AC_XLD_D2_RK1(rng, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_DM_AC_XLD_D2_RK1
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
    impure module function getMNR_RNGX_AM_AC_UXD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D2_RK5
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
    impure module function getMNR_RNGX_AM_AC_UXD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D2_RK4
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
    impure module function getMNR_RNGX_AM_AC_UXD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D2_RK3
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
    impure module function getMNR_RNGX_AM_AC_UXD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D2_RK2
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
    impure module function getMNR_RNGX_AM_AC_UXD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_UXD_D2_RK1
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
    impure module function getMNR_RNGX_AM_AC_XLD_D2_RK5(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D2_RK5
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
    impure module function getMNR_RNGX_AM_AC_XLD_D2_RK4(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D2_RK4
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
    impure module function getMNR_RNGX_AM_AC_XLD_D2_RK3(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D2_RK3
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
    impure module function getMNR_RNGX_AM_AC_XLD_D2_RK2(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D2_RK2
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
    impure module function getMNR_RNGX_AM_AC_XLD_D2_RK1(rng, mean, chol, subset, nsam) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMNR_RNGX_AM_AC_XLD_D2_RK1
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
    !>  Return a (collection) of random vector(s) of size `ndim` from the `ndim`-dimensional MultiVariate Normal (MVN) distribution,
    !>  optionally with the specified input `mean(1:ndim)` and the specified `subset` of the Cholesky Factorization of the Covariance matrix of the MVN distribution.
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
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              containing the Multivariate Normal distributed random output vector.<br>
    !>  \param[in]      mean    :   The input `contiguous` vector of shape `(1:ndim)`, of the same type and kind as the output `rand`, representing the mean of the Multivariate Normal distribution.<br>
    !>                              (**optional**, default = `[(0., i = 1, size(rand))]`. It must be present if the input argument `chol` is missing.)
    !>  \param[in]      chol    :   The input `contiguous` matrix of shape `(ndim, ndim)` whose specified triangular `subset` contains the [Cholesky Factorization](@ref pm_matrixChol) of the covariance matrix of the MVN distribution.<br>
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
    !>  \interface{setMultiNormRand}
    !>  \code{.F90}
    !>
    !>      use pm_distMultiNorm, only: setMultiNormRand
    !>
    !>      ! single vector, using default rng
    !>
    !>      call setMultiNormRand(rand(1:ndim), mean(1:ndim))
    !>      call setMultiNormRand(rand(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>      call setMultiNormRand(rand(1:ndim), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! single vector, using custom rng
    !>
    !>      call setMultiNormRand(rng, rand(1:ndim), mean(1:ndim))
    !>      call setMultiNormRand(rng, rand(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>      call setMultiNormRand(rng, rand(1:ndim), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using default rng
    !>
    !>      call setMultiNormRand(rand(1:ndim, 1:nsam), mean(1:ndim))
    !>      call setMultiNormRand(rand(1:ndim, 1:nsam), chol(1:ndim, 1:ndim), subset)
    !>      call setMultiNormRand(rand(1:ndim, 1:nsam), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
    !>
    !>      ! collection of `nsam` vectors, using custom rng
    !>
    !>      call setMultiNormRand(rng, rand(1:ndim, 1:nsam), mean(1:ndim))
    !>      call setMultiNormRand(rng, rand(1:ndim, 1:nsam), chol(1:ndim, 1:ndim), subset)
    !>      call setMultiNormRand(rng, rand(1:ndim, 1:nsam), mean(1:ndim), chol(1:ndim, 1:ndim), subset)
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
    !>  \see
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getMultiNormLogPDF](@ref pm_distMultiNorm::getMultiNormLogPDF)<br>
    !>
    !>  \example{setMultiNormRand}
    !>  \include{lineno} example/pm_distMultiNorm/setMultiNormRand/main.F90
    !>  \compilef{setMultiNormRand}
    !>  \output{setMultiNormRand}
    !>  \include{lineno} example/pm_distMultiNorm/setMultiNormRand/main.out.F90
    !>  \postproc{setMultiNormRand}
    !>  \include{lineno} example/pm_distMultiNorm/setMultiNormRand/main.py
    !>  \vis{setMultiNormRand}
    !>  \image html pm_distMultiNorm/setMultiNormRand/setMultiNormRand.RK.png width=700
    !>  \image html pm_distMultiNorm/setMultiNormRand/setMultiNormRandMean.RK.png width=700
    !>  \image html pm_distMultiNorm/setMultiNormRand/setMultiNormRandChol.RK.png width=700
    !>  \image html pm_distMultiNorm/setMultiNormRand/setMultiNormRandMeanChol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distMultiNorm](@ref test_pm_distMultiNorm)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setMNR_RNGD_DM_DC_XXX_D1_RK5()
    !>             |||| || || ||| || |||
    !>             |||| || || ||| || |||
    !>             |||| || || ||| || |||
    !>             |||| || || ||| || |The Kind of the output array.
    !>             |||| || || ||| || The type of the output array: R => Real.
    !>             |||| || || ||| The Dimension of the output array.
    !>             |||| || || The subset of the Cholesky factor: D/U/L => Default/Upper/Lower.
    !>             |||| || The Cholesky factor of the covariance matrix of the distribution: DC/AC => Default/Arbitrary Cholesky
    !>             |||| The mean of the distribution: DM/AM => Default/Arbitrary Mean
    !>             The random number generator: RNG D/F/X => Default/Fortran/Xoroshiro256++
    !>  \endcode
    !>
    !>  \todo
    !>  \pmed
    !>  The access pattern for the upper-diagonal subset of `chol` is non `contiguous` in the current implementation.<br>
    !>  The access pattern can be likely made `contiguous` by an appropriate implementation.<br>
    !>
    !>  \final{setMultiNormRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 12:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

    ! D1 RNGD

    interface setMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGD_DM_DC_XXX_D1_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGD_DM_DC_XXX_D1_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGD_DM_DC_XXX_D1_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGD_DM_DC_XXX_D1_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGD_DM_DC_XXX_D1_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGD_AM_DC_XXX_D1_RK5(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGD_AM_DC_XXX_D1_RK4(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGD_AM_DC_XXX_D1_RK3(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGD_AM_DC_XXX_D1_RK2(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGD_AM_DC_XXX_D1_RK1(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_UXD_D1_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_UXD_D1_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_UXD_D1_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_UXD_D1_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_UXD_D1_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_XLD_D1_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_XLD_D1_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_XLD_D1_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_XLD_D1_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGD_DM_AC_XLD_D1_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_UXD_D1_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_UXD_D1_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_UXD_D1_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_UXD_D1_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_UXD_D1_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_XLD_D1_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_XLD_D1_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_XLD_D1_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_XLD_D1_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGD_AM_AC_XLD_D1_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D1_RK1
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

    interface setMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGF_DM_DC_XXX_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGF_DM_DC_XXX_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGF_DM_DC_XXX_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGF_DM_DC_XXX_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGF_DM_DC_XXX_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGF_AM_DC_XXX_D1_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGF_AM_DC_XXX_D1_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGF_AM_DC_XXX_D1_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGF_AM_DC_XXX_D1_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGF_AM_DC_XXX_D1_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_UXD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_UXD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_UXD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_UXD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_UXD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D1_RK1
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
    impure module subroutine setMNR_RNGF_DM_AC_XLD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_XLD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_XLD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_XLD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGF_DM_AC_XLD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D1_RK1
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
    impure module subroutine setMNR_RNGF_AM_AC_UXD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_UXD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_UXD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_UXD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_UXD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D1_RK1
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
    impure module subroutine setMNR_RNGF_AM_AC_XLD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_XLD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_XLD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_XLD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setMNR_RNGF_AM_AC_XLD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D1_RK1
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

    interface setMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D1_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D1_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D1_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D1_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D1_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D1_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D1_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D1_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D1_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D1_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D1_RK1
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
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D1_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D1_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D1_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D1_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D1_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D1_RK1
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
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D1_RK1
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
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D1_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D1_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D1_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D1_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D1_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D1_RK1
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

    interface setMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGD_DM_DC_XXX_D2_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGD_DM_DC_XXX_D2_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGD_DM_DC_XXX_D2_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGD_DM_DC_XXX_D2_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGD_DM_DC_XXX_D2_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGD_AM_DC_XXX_D2_RK5(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGD_AM_DC_XXX_D2_RK4(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGD_AM_DC_XXX_D2_RK3(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGD_AM_DC_XXX_D2_RK2(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGD_AM_DC_XXX_D2_RK1(rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGD_DM_AC_UXD_D2_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGD_DM_AC_UXD_D2_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGD_DM_AC_UXD_D2_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGD_DM_AC_UXD_D2_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGD_DM_AC_UXD_D2_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGD_DM_AC_XLD_D2_RK5(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGD_DM_AC_XLD_D2_RK4(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGD_DM_AC_XLD_D2_RK3(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGD_DM_AC_XLD_D2_RK2(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGD_DM_AC_XLD_D2_RK1(rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_DM_AC_XLD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGD_AM_AC_UXD_D2_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGD_AM_AC_UXD_D2_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGD_AM_AC_UXD_D2_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGD_AM_AC_UXD_D2_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGD_AM_AC_UXD_D2_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_UXD_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGD_AM_AC_XLD_D2_RK5(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGD_AM_AC_XLD_D2_RK4(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGD_AM_AC_XLD_D2_RK3(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGD_AM_AC_XLD_D2_RK2(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGD_AM_AC_XLD_D2_RK1(rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGD_AM_AC_XLD_D2_RK1
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

    interface setMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGF_DM_DC_XXX_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGF_DM_DC_XXX_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGF_DM_DC_XXX_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGF_DM_DC_XXX_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGF_DM_DC_XXX_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGF_AM_DC_XXX_D2_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGF_AM_DC_XXX_D2_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGF_AM_DC_XXX_D2_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGF_AM_DC_XXX_D2_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGF_AM_DC_XXX_D2_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMNR_RNGF_DM_AC_UXD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGF_DM_AC_UXD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGF_DM_AC_UXD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGF_DM_AC_UXD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGF_DM_AC_UXD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_UXD_D2_RK1
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
    module subroutine setMNR_RNGF_DM_AC_XLD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGF_DM_AC_XLD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGF_DM_AC_XLD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGF_DM_AC_XLD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGF_DM_AC_XLD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_DM_AC_XLD_D2_RK1
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
    module subroutine setMNR_RNGF_AM_AC_UXD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGF_AM_AC_UXD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGF_AM_AC_UXD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGF_AM_AC_UXD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGF_AM_AC_UXD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_UXD_D2_RK1
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
    module subroutine setMNR_RNGF_AM_AC_XLD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMNR_RNGF_AM_AC_XLD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMNR_RNGF_AM_AC_XLD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMNR_RNGF_AM_AC_XLD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMNR_RNGF_AM_AC_XLD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGF_AM_AC_XLD_D2_RK1
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

    interface setMultiNormRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D2_RK5(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D2_RK4(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D2_RK3(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D2_RK2(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_DM_DC_XXX_D2_RK1(rng, rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D2_RK5(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D2_RK4(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D2_RK3(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D2_RK2(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_AM_DC_XXX_D2_RK1(rng, rand, mean)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_DC_XXX_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_UXD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_UXD_D2_RK1
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
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D2_RK5(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D2_RK4(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D2_RK3(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D2_RK2(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_DM_AC_XLD_D2_RK1(rng, rand, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_DM_AC_XLD_D2_RK1
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
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_UXD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_UXD_D2_RK1
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
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D2_RK5(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D2_RK4(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D2_RK3(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D2_RK2(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKG)               , intent(out)   , contiguous    :: rand(:,:)
        real(RKG)               , intent(in)    , contiguous    :: mean(:), chol(:,:)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMNR_RNGX_AM_AC_XLD_D2_RK1(rng, rand, mean, chol, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMNR_RNGX_AM_AC_XLD_D2_RK1
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

end module pm_distMultiNorm