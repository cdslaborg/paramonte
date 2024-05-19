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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Lognormal distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Lognormal distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the Lognormal distribution is defined with the two location and scale parameters \f$(\mu \in (-\infty, +\infty), \sigma > 0)\f$ as,
!>  \f{equation}{
!>      \large
!>      \pi(x | \mu, \sigma) =
!>      \frac{1}{x\sigma\sqrt{2\pi}}\exp\bigg( -\frac{\big(\log(x) - \mu\big)^2}{2\sigma^2} \bigg) ~,~ x \in (0, +\infty) ~.
!>  \f}
!>
!>  The **CDF** of the Lognormal distribution is defined with the two location and scale parameters \f$(\mu \in (-\infty, +\infty), \sigma > 0)\f$ as,
!>  \f{equation}{
!>      \large
!>      \mathrm{CDF}(x | \mu, \sigma) =
!>      \frac{1}{2} \bigg[ 1 + \mathrm{erf} \bigg( \frac{\log(x) - \mu}{\sigma\sqrt{2}} \bigg) \bigg] ~,~ x \in (0, +\infty) ~.
!>  \f}
!>
!>  \test
!>  [test_pm_distLogNorm](@ref test_pm_distLogNorm)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distLogNorm

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter         :: MODULE_NAME = "@pm_distLogNorm"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Lognormal
    !>  as defined in the description of [pm_distLogNorm](@ref pm_distLogNorm).
    !>
    !>  \details
    !>  See the documentation of [pm_distLogNorm](@ref pm_distLogNorm) for the definition of the Lognormal distribution.
    !>
    !>  \interface{distLogNorm_type}
    !>  \code{.F90}
    !>
    !>      use pm_distLogNorm, only: distLogNorm_type
    !>      type(distLogNorm_type) :: distLogNorm
    !>
    !>      distLogNorm = distLogNorm_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distLogNorm](@ref test_pm_distLogNorm)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distLogNorm_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distLogNorm_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the natural logarithm of probability density function (PDF) of the univariate Lognormal distribution.
    !>
    !>  \param[in]  x           :   The input positive-valued scalar or array of the same shape as other array-like arguments,
    !>                              of type `real` of kind \RKALL, representing the point(s) at which the PDF must be computed.
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as `x` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  sigma       :   The input scalar of the same type and kind as `x` representing the inverse of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar or array of the same shape as the input array-like arguments, of the same
    !>                              type and kind as the input `x` representing the PDF of the Lognormal distribution at `x`.
    !>
    !>  \interface{getLogNormLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distLogNorm, only: getLogNormLogPDF
    !>
    !>      logPDF = getLogNormLogPDF(x, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `x > 0.` and `0. < sigma` must hold for the corresponding procedure arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  See [setLogNormLogPDF](@ref pm_distLogNorm::setLogNormLogPDF) for a more performant but less flexible interface of the same functionality.<br>
    !>
    !>  \see
    !>  [setLogNormLogPDF](@ref pm_distLogNorm::setLogNormLogPDF)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [setNormLogPDF](@ref pm_distNorm::setNormLogPDF)<br>
    !>
    !>  \example{getLogNormLogPDF}
    !>  \include{lineno} example/pm_distLogNorm/getLogNormLogPDF/main.F90
    !>  \compilef{getLogNormLogPDF}
    !>  \output{getLogNormLogPDF}
    !>  \include{lineno} example/pm_distLogNorm/getLogNormLogPDF/main.out.F90
    !>  \postproc{getLogNormLogPDF}
    !>  \include{lineno} example/pm_distLogNorm/getLogNormLogPDF/main.py
    !>  \vis{getLogNormLogPDF}
    !>  \image html pm_distLogNorm/getLogNormLogPDF/getLogNormLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogNorm](@ref test_pm_distLogNorm)
    !>
    !>  \final{getLogNormLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogNormLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogNormLogPDF_RK5(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormLogPDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogNormLogPDF_RK4(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormLogPDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogNormLogPDF_RK3(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormLogPDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogNormLogPDF_RK2(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormLogPDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogNormLogPDF_RK1(x, mu, sigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormLogPDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(in)                    :: x
        real(RKG)              , intent(in) , optional         :: mu, sigma
        real(RKG)                                              :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLogNormLogPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the natural logarithm of probability density function (PDF) of the univariate Lognormal distribution.
    !>
    !>  \param[out] logPDF      :   The output scalar or array of the same shape as the input array-like arguments, of the same type
    !>                              and kind as the input `logx` representing the PDF of the Lognormal distribution at `exp(logx)`.
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of type `real` of kind \RKALL, representing the natural logarithms of the point(s) at which the PDF must be computed (i.e., \f$\log(x)\f$).
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `logx` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  invSigma    :   The input scalar of the same type and kind as `logx` representing the inverse of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`, must be present <b>if and only if</b> `logInvSigma` is also present)
    !>  \param[in]  logInvSigma :   The input scalar of the same type and kind as `logx` representing the natural logarithm of the inverse of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `0`, must be present <b>if and only if</b> `invSigma` is also present).
    !>
    !>  \interface{setLogNormLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distLogNorm, only: setLogNormLogPDF
    !>
    !>      call setLogNormLogPDF(logPDF, logx)
    !>      call setLogNormLogPDF(logPDF, logx, mu)
    !>      call setLogNormLogPDF(logPDF, logx, invSigma, logInvSigma)
    !>      call setLogNormLogPDF(logPDF, logx, mu, invSigma, logInvSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `invSigma > 0.` must hold for the corresponding procedure argument.<br>
    !>  The condition `log(invSigma) == logInvSigma` must hold within a small range of computer precision for the corresponding procedure arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  See [getLogNormLogPDF](@ref pm_distLogNorm::getLogNormLogPDF) for a less performant but more flexible interface of the same functionality.<br>
    !>
    !>  \see
    !>  [getLogNormLogPDF](@ref pm_distLogNorm::getLogNormLogPDF)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [setNormLogPDF](@ref pm_distNorm::setNormLogPDF)<br>
    !>
    !>  \example{setLogNormLogPDF}
    !>  \include{lineno} example/pm_distLogNorm/setLogNormLogPDF/main.F90
    !>  \compilef{setLogNormLogPDF}
    !>  \output{setLogNormLogPDF}
    !>  \include{lineno} example/pm_distLogNorm/setLogNormLogPDF/main.out.F90
    !>  \postproc{setLogNormLogPDF}
    !>  \include{lineno} example/pm_distLogNorm/setLogNormLogPDF/main.py
    !>  \vis{setLogNormLogPDF}
    !>  \image html pm_distLogNorm/setLogNormLogPDF/setLogNormLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogNorm](@ref test_pm_distLogNorm)
    !>
    !>  \todo
    !>  \pmed
    !>  A performant vectorized `logPDF(:)` version of the subroutines under this generic interface could be added in the future.
    !>
    !>  \final{setLogNormLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setLogNormLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDD_RK5(logPDF, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDD_RK4(logPDF, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDD_RK3(logPDF, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDD_RK2(logPDF, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDD_RK1(logPDF, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMD_RK5(logPDF, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMD_RK4(logPDF, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMD_RK3(logPDF, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMD_RK2(logPDF, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMD_RK1(logPDF, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDS_RK5(logPDF, logx, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDS_RK4(logPDF, logx, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDS_RK3(logPDF, logx, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDS_RK2(logPDF, logx, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormLogPDFDS_RK1(logPDF, logx, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFDS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: invSigma, logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMS_RK5(logPDF, logx, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMS_RK4(logPDF, logx, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMS_RK3(logPDF, logx, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMS_RK2(logPDF, logx, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormLogPDFMS_RK1(logPDF, logx, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormLogPDFMS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)              , intent(out)                   :: logPDF
        real(RKG)              , intent(in)                    :: logx
        real(RKG)              , intent(in)                    :: mu
        real(RKG)              , intent(in)                    :: invSigma
        real(RKG)              , intent(in)                    :: logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setLogNormLogPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the univariate Lognormal distribution.
    !>
    !>  \param[in]  x           :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of type `real` of kind \RKALL, representing the point(s) at which the CDF must be computed.
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as `x` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  sigma       :   The input scalar of the same type and kind as `x` representing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as the input array-like arguments, of the same type
    !>                              and kind as the input `x` representing the CDF of the Lognormal distribution at the given input `x`.
    !>
    !>  \interface{getLogNormCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distLogNorm, only: getLogNormCDF
    !>
    !>      cdf = getLogNormCDF(x, mu = mu, sigma = sigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < sigma` must hold for the corresponding procedure argument.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \example{getLogNormCDF}
    !>  \include{lineno} example/pm_distLogNorm/getLogNormCDF/main.F90
    !>  \compilef{getLogNormCDF}
    !>  \output{getLogNormCDF}
    !>  \include{lineno} example/pm_distLogNorm/getLogNormCDF/main.out.F90
    !>  \postproc{getLogNormCDF}
    !>  \include{lineno} example/pm_distLogNorm/getLogNormCDF/main.py
    !>  \vis{getLogNormCDF}
    !>  \image html pm_distLogNorm/getLogNormCDF/getLogNormCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogNorm](@ref test_pm_distLogNorm)
    !>
    !>  \final{getLogNormCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogNormCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogNormCDF_RK5(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogNormCDF_RK4(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogNormCDF_RK3(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogNormCDF_RK2(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogNormCDF_RK1(x, mu, sigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogNormCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)                            :: x
        real(RKG)           , intent(in)    , optional              :: mu
        real(RKG)           , intent(in)    , optional              :: sigma
        real(RKG)                                                   :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLogNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the univariate Lognormal distribution.
    !>
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of type `real` of kind \RKALL, representing the natural logarithm of the point(s) at which the CDF must be computed.
    !>  \param[in]  mu          :   The input scalar or array of the same shape as other array-like arguments of the same type
    !>                              and kind as `logx` representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  invSigma    :   The input scalar of the same type and kind as `logx` representing the inverse of the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`. It can be present <b>if and only if</b> is also present.)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as the input array-like arguments, of the same type
    !>                              and kind as the input `logx` representing the CDF of the Lognormal distribution at the specified point `exp(logx)`.
    !>
    !>  \interface{setLogNormCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distLogNorm, only: setLogNormCDF
    !>
    !>      call setLogNormCDF(cdf, logx)
    !>      call setLogNormCDF(cdf, logx, mu)
    !>      call setLogNormCDF(cdf, logx, mu, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `invSigma > 0.` must hold for the corresponding procedure argument.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \example{setLogNormCDF}
    !>  \include{lineno} example/pm_distLogNorm/setLogNormCDF/main.F90
    !>  \compilef{setLogNormCDF}
    !>  \output{setLogNormCDF}
    !>  \include{lineno} example/pm_distLogNorm/setLogNormCDF/main.out.F90
    !>  \postproc{setLogNormCDF}
    !>  \include{lineno} example/pm_distLogNorm/setLogNormCDF/main.py
    !>  \vis{setLogNormCDF}
    !>  \image html pm_distLogNorm/setLogNormCDF/setLogNormCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogNorm](@ref test_pm_distLogNorm)
    !>
    !>  \final{setLogNormCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setLogNormCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormCDFDD_RK5(cdf, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormCDFDD_RK4(cdf, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormCDFDD_RK3(cdf, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormCDFDD_RK2(cdf, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(out)                           :: cdf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormCDFDD_RK1(cdf, logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormCDFMD_RK5(cdf, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormCDFMD_RK4(cdf, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormCDFMD_RK3(cdf, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormCDFMD_RK2(cdf, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormCDFMD_RK1(cdf, logx, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogNormCDFMS_RK5(cdf, logx, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMS_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogNormCDFMS_RK4(cdf, logx, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMS_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogNormCDFMS_RK3(cdf, logx, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMS_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogNormCDFMS_RK2(cdf, logx, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMS_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogNormCDFMS_RK1(cdf, logx, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogNormCDFMS_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)                           :: cdf
        real(RKG)           , intent(in)                            :: logx
        real(RKG)           , intent(in)                            :: mu
        real(RKG)           , intent(in)                            :: invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setLogNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distLogNorm