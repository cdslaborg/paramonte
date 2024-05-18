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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Exponential distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Exponential distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the **Exponential distribution** with the two location and scale parameters \f$(\mu, \sigma)\f$ is defined as,<br>
!>
!>  \f{equation}{
!>      \large
!>      F(x | \mu, \sigma) =
!>      \begin{cases}
!>          \frac{1}{\sigma} \exp(-\frac{x - \mu}{\sigma})  &,~ \mu \leq x < +\infty \\
!>          0                                               &,~ x < \mu
!>      \end{cases}
!>  \f}
!>
!>  where \f$-\infty < \mu < +\infty\f$ is the location parameter and \f$0 < \sigma < +\infty\f$ is the scale parameter of the distribution.<br>
!>
!>  The corresponding **CDF** with the two location and scale parameters \f$(\mu, \sigma)\f$ is defined as,<br>
!>  \f{equation}{
!>      \large
!>      \mathrm{CDF}(x | \mu, \sigma) =
!>      \begin{cases}
!>          1 - \exp(-\frac{x - \mu}{\sigma})   &,~ \mu \leq x < +\infty \\
!>          0                                   &,~ x < \mu
!>      \end{cases}
!>  \f}
!>  where \f$-\infty < \mu < +\infty\f$ is the location parameter and \f$0 < \sigma < +\infty\f$ is the scale parameter of the distribution.<br>
!>
!>  <b>Random Number Generation</b>
!>
!>  Assuming \f$U \in (0, 1]\f$ is a uniformly-distributed random variate,<br>
!>  \f{equation}{
!>      \large
!>      T = -\sigma\log(U) + \mu ~,
!>  \f}
!>  is a random variate from the Exponential distribution with the location parameter \f$-\infty < \mu < +\infty\f$
!>  and the scale parameter \f$0 < \sigma < +\infty\f$ such that \f$\mu \leq T < +\infty\f$.<br>
!>
!>  \note
!>  Note that the mean of the Exponential distribution is the scale parameter: \f$\mathrm{mean} = \sigma\f$).<br>
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>  [pm_distNegExp](@ref pm_distNegExp)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getExpLogPDF_vs_setExpLogPDF, The runtime performance of [getExpLogPDF](@ref pm_distExp::getExpLogPDF) vs. [setExpLogPDF](@ref pm_distExp::setExpLogPDF)}
!>  \include{lineno} benchmark/pm_distExp/getExpLogPDF_vs_setExpLogPDF/main.F90
!>  \compilefb{getExpLogPDF_vs_setExpLogPDF}
!>  \postprocb{getExpLogPDF_vs_setExpLogPDF}
!>  \include{lineno} benchmark/pm_distExp/getExpLogPDF_vs_setExpLogPDF/main.py
!>  \visb{getExpLogPDF_vs_setExpLogPDF}
!>  \image html benchmark/pm_distExp/getExpLogPDF_vs_setExpLogPDF/benchmark.getExpLogPDF_vs_setExpLogPDF.runtime.png width=1000
!>  \image html benchmark/pm_distExp/getExpLogPDF_vs_setExpLogPDF/benchmark.getExpLogPDF_vs_setExpLogPDF.runtime.ratio.png width=1000
!>  \moralb{getExpLogPDF_vs_setExpLogPDF}
!>      -#  The procedures under the generic interface [getExpLogPDF](@ref pm_distExp::getExpLogPDF) are functions while
!>          the procedures under the generic interface [setExpLogPDF](@ref pm_distExp::setExpLogPDF) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than
!>          the subroutine interface when the input `array` size is small.<br>
!>          Otherwise, the difference appears to be marginal and insignificant in most practical situations.<br>
!>
!>  \benchmark{setExpLogPDF-logInvSigma-missing_vs_present, The runtime performance of [setExpLogPDF](@ref pm_distExp::setExpLogPDF) with and without `logInvSigma`}
!>  \include{lineno} benchmark/pm_distExp/setExpLogPDF-logInvSigma-missing_vs_present/main.F90
!>  \compilefb{setExpLogPDF-logInvSigma-missing_vs_present}
!>  \postprocb{setExpLogPDF-logInvSigma-missing_vs_present}
!>  \include{lineno} benchmark/pm_distExp/setExpLogPDF-logInvSigma-missing_vs_present/main.py
!>  \visb{setExpLogPDF-logInvSigma-missing_vs_present}
!>  \image html benchmark/pm_distExp/setExpLogPDF-logInvSigma-missing_vs_present/benchmark.setExpLogPDF-logInvSigma-missing_vs_present.runtime.png width=1000
!>  \image html benchmark/pm_distExp/setExpLogPDF-logInvSigma-missing_vs_present/benchmark.setExpLogPDF-logInvSigma-missing_vs_present.runtime.ratio.png width=1000
!>  \moralb{setExpLogPDF-logInvSigma-missing_vs_present}
!>      -#  The procedures under the generic interface [setExpLogPDF](@ref pm_distExp::setExpLogPDF)
!>          accept an extra argument `logInvSigma = log(invSigma)` while the procedures under the generic interface
!>          [getExpLogPDF](@ref pm_distExp::getExpLogPDF) compute this term internally with every procedure call.<br>
!>          In the presence of this argument, the logarithmic computation `log(invSigma)` will be avoided.<br>
!>          As such, the presence of `logInvSigma` is expected to lead to faster computations.<br>
!>
!>  \test
!>  [test_pm_distExp](@ref test_pm_distExp)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distExp

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Exponential
    !>  as defined in the description of [pm_distExp](@ref pm_distExp).
    !>
    !>  \details
    !>  See the documentation of [pm_distExp](@ref pm_distExp) for the definition of the Exponential distribution.
    !>
    !>  \interface{distExp_type}
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: distExp_type
    !>      type(distExp_type) :: distExp
    !>
    !>      distExp = distExp_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distExp_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distExp_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the
    !>  <b>Exponential distribution</b> for an input `x` within the support of the distribution \f$[\mu, +\infty)\f$.
    !>
    !>  \param[in]  x           :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the values at which the PDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  invSigma    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the inverse scale (i.e., the rate or the inverse mean) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `logPDF`                :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`,
    !>                              containing the natural logarithm of the PDF of the distribution at the specified point.<br>
    !>
    !>  \interface{getExpLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: getExpLogPDF
    !>
    !>      logPDF = getExpLogPDF(x, mu = mu, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `mu <= x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setExpLogPDF](@ref pm_distExp::setExpLogPDF)<br>
    !>
    !>  \example{getExpLogPDF}
    !>  \include{lineno} example/pm_distExp/getExpLogPDF/main.F90
    !>  \compilef{getExpLogPDF}
    !>  \output{getExpLogPDF}
    !>  \include{lineno} example/pm_distExp/getExpLogPDF/main.out.F90
    !>  \postproc{getExpLogPDF}
    !>  \include{lineno} example/pm_distExp/getExpLogPDF/main.py
    !>  \vis{getExpLogPDF}
    !>  \image html pm_distExp/getExpLogPDF/getExpLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getExpLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getExpLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpLogPDFXMI_RK5(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpLogPDFXMI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpLogPDFXMI_RK4(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpLogPDFXMI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpLogPDFXMI_RK3(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpLogPDFXMI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpLogPDFXMI_RK2(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpLogPDFXMI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpLogPDFXMI_RK1(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpLogPDFXMI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the
    !>  <b>Exponential distribution</b> for an input `x` within the support of the distribution \f$[\mu, +\infty)\f$.
    !>
    !>  \param[out] logPDF      :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`,
    !>                              containing the natural logarithm of the PDF of the distribution at the specified point.<br>
    !>  \param[in]  x           :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the values at which the PDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  invSigma    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the inverse scale (i.e., the rate or the inverse mean) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`. It must be present <b>if and only if</b> `logInvSigma` is also present.)
    !>  \param[in]  logInvSigma :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the natural logarithm of the inverse scale parameter `invSigma`.<br>
    !>                              If present, it will lead to significantly faster computations of the PDF.<br>
    !>                              (**optional**, default = `log(invSigma)`. It must be present <b>if and only if</b> `invSigma` is also present.)
    !>
    !>  \interface{setExpLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: setExpLogPDF
    !>
    !>      call setExpLogPDF(logPDF, x)
    !>      call setExpLogPDF(logPDF, x, mu)
    !>      call setExpLogPDF(logPDF, x, invSigma, logInvSigma)
    !>      call setExpLogPDF(logPDF, x, mu, invSigma, logInvSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `mu <= x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  The condition `logInvSigma = log(invSigma)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getExpLogPDF](@ref pm_distExp::getExpLogPDF)<br>
    !>
    !>  \example{setExpLogPDF}
    !>  \include{lineno} example/pm_distExp/setExpLogPDF/main.F90
    !>  \compilef{setExpLogPDF}
    !>  \output{setExpLogPDF}
    !>  \include{lineno} example/pm_distExp/setExpLogPDF/main.out.F90
    !>  \postproc{setExpLogPDF}
    !>  \include{lineno} example/pm_distExp/setExpLogPDF/main.py
    !>  \vis{setExpLogPDF}
    !>  \image html pm_distExp/setExpLogPDF/setExpLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setExpLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setExpLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpLogPDFXDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpLogPDFXDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpLogPDFXDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpLogPDFXDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpLogPDFXDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpLogPDFXMD_RK5(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpLogPDFXMD_RK4(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpLogPDFXMD_RK3(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpLogPDFXMD_RK2(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpLogPDFXMD_RK1(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpLogPDFXDI_RK5(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpLogPDFXDI_RK4(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpLogPDFXDI_RK3(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpLogPDFXDI_RK2(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpLogPDFXDI_RK1(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXDI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpLogPDFXMI_RK5(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpLogPDFXMI_RK4(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpLogPDFXMI_RK3(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpLogPDFXMI_RK2(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpLogPDFXMI_RK1(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpLogPDFXMI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logPDF
        real(RKC)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  <b>Exponential distribution</b> for an input `x` within the support of the distribution \f$[\mu, +\infty)\f$.
    !>
    !>  \param[in]  x           :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the values at which the CDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`)
    !>  \param[in]  invSigma    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the inverse scale (i.e., the rate or the inverse mean) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`)
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`,
    !>                              containing the CDF of the distribution at the specified point.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: getExpCDF
    !>
    !>      cdf = getExpCDF(x, mu = mu, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `mu <= x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setExpCDF](@ref pm_distExp::setExpCDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExp/getExpCDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExp/getExpCDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExp/getExpCDF/main.py
    !>  \vis
    !>  \image html pm_distExp/getExpCDF/getExpCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getExpCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getExpCDFXMI_RK5(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpCDFXMI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getExpCDFXMI_RK4(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpCDFXMI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getExpCDFXMI_RK3(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpCDFXMI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getExpCDFXMI_RK2(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpCDFXMI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getExpCDFXMI_RK1(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpCDFXMI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x
        real(RKC)   , intent(in)    , optional  :: mu, invSigma
        real(RKC)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  <b>Exponential distribution</b> for an input `x` within the support of the distribution \f$[\mu, +\infty)\f$.
    !>
    !>  \param[out] cdf         :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`,
    !>                              containing the CDF of the distribution at the specified point.<br>
    !>  \param[in]  x           :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of <br>
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              containing the values at which the CDF must be computed.<br>
    !>  \param[in]  mu          :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`. If present, then `invSigma` must also be present.)
    !>  \param[in]  invSigma    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the inverse scale (i.e., the rate or the inverse mean) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`.)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: setExpCDF
    !>
    !>      call setExpCDF(cdf, x)
    !>      call setExpCDF(cdf, x, invSigma)
    !>      call setExpCDF(cdf, x, mu, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `mu <= x` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getExpCDF](@ref pm_distExp::getExpCDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distExp/setExpCDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distExp/setExpCDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distExp/setExpCDF/main.py
    !>  \vis
    !>  \image html pm_distExp/setExpCDF/setExpCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setExpCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpCDFXDD_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpCDFXDD_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpCDFXDD_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpCDFXDD_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpCDFXDD_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpCDFXDI_RK5(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpCDFXDI_RK4(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpCDFXDI_RK3(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpCDFXDI_RK2(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpCDFXDI_RK1(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXDI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpCDFXMI_RK5(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXMI_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpCDFXMI_RK4(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXMI_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpCDFXMI_RK3(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXMI_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpCDFXMI_RK2(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXMI_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpCDFXMI_RK1(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpCDFXMI_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Exponential distribution,
    !>  optionally with the specified input location and scale parameters `mu, sigma`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distExp](@ref pm_distExp) for more details.<br>
    !>
    !>  \param[in]  sigma   :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                          <ul>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          representing the scale parameter of the distribution.<br>
    !>  \param[in]  mu      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `sigma`,
    !>                          representing the location parameter of the distribution.<br>
    !>                          (**optional**, default = `0`.)
    !>
    !>  \return
    !>  `rand`              :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `sigma`,
    !>                          containing random value(s) from the distribution.<br>
    !>
    !>  \interface{getExpRand}
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: getExpRand
    !>
    !>      rand = getExpRand(sigma, mu = mu)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < sigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setExpRand](@ref pm_distExp::setExpRand)<br>
    !>
    !>  \example{getExpRand}
    !>  \include{lineno} example/pm_distExp/getExpRand/main.F90
    !>  \compilef{getExpRand}
    !>  \output{getExpRand}
    !>  \include{lineno} example/pm_distExp/getExpRand/main.out.F90
    !>  \postproc{getExpRand}
    !>  \include{lineno} example/pm_distExp/getExpRand/main.py
    !>  \vis{getExpRand}
    !>  \image html pm_distExp/getExpRand/getExpRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \final{getExpRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getExpRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getExpRandSM_RK5(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpRandSM_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                                   :: rand
        real(RKC)   , intent(in)                    :: sigma
        real(RKC)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getExpRandSM_RK4(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpRandSM_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                                   :: rand
        real(RKC)   , intent(in)                    :: sigma
        real(RKC)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getExpRandSM_RK3(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpRandSM_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                                   :: rand
        real(RKC)   , intent(in)                    :: sigma
        real(RKC)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getExpRandSM_RK2(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpRandSM_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                                   :: rand
        real(RKC)   , intent(in)                    :: sigma
        real(RKC)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getExpRandSM_RK1(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExpRandSM_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                                   :: rand
        real(RKC)   , intent(in)                    :: sigma
        real(RKC)   , intent(in)    , optional      :: mu
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Exponential distribution,
    !>  optionally with the specified input location and scale parameters `mu, sigma`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distExp](@ref pm_distExp) for more details.<br>
    !>
    !>  \param[inout]   rand    :   The input/output scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ol>
    !>                              that,
    !>                              <ol>
    !>                                  <li>    must contain uniformly-distributed random value(s) in range \f$[0, 1)\f$ on input, and
    !>                                  <li>    will contain exponentially-distributed random value(s) in range \f$[0, +\infty)\f$ on output.
    !>                              </ol>
    !>                              The uniformly-distributed random value(s) can be readily obtained via
    !>                              <ol>
    !>                                  <li>    the Fortran intrinsic procedure `random_number()` or,
    !>                                  <li>    via [getUnifRand()](@ref pm_distUnif::getUnifRand).
    !>                              </ol>
    !>                              Supplying this argument as an input ensures the purity of the procedures, allowing further compiler optimizations.<br>
    !>  \param[in]      sigma   :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `rand`,
    !>                              representing the scale parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`.)
    !>  \param[in]      mu      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `rand`,
    !>                              representing the location parameter of the distribution.<br>
    !>                              (**optional**, default = `0`. It can be present <b>if and only if</b> `sigma` is also present.)
    !>
    !>  \interface{setExpRand}
    !>  \code{.F90}
    !>
    !>      use pm_distExp, only: setExpRand
    !>
    !>      call setExpRand(rand)
    !>      call setExpRand(rand, sigma)
    !>      call setExpRand(rand, sigma, mu)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= rand < 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getExpRand](@ref pm_distExp::getExpRand)<br>
    !>
    !>  \example{setExpRand}
    !>  \include{lineno} example/pm_distExp/setExpRand/main.F90
    !>  \compilef{setExpRand}
    !>  \output{setExpRand}
    !>  \include{lineno} example/pm_distExp/setExpRand/main.out.F90
    !>  \postproc{setExpRand}
    !>  \include{lineno} example/pm_distExp/setExpRand/main.py
    !>  \vis{setExpRand}
    !>  \image html pm_distExp/setExpRand/setExpRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distExp](@ref test_pm_distExp)
    !>
    !>  \todo
    !>  \pvhigh
    !>  The Intel ifort bug (described below) appears to have been resolved in ifort 2021.4.<br>
    !>  Therefore, once TACC and other relevant supercomputers have ifort >=2021.4 installed, the bug workaround must be resolved.<br>
    !>
    !>  \bug
    !>  \status possibly \resolved in \ifort{2021.4}
    !>  \source \ifort{2021.2}
    !>  \desc
    !>  The \ifort{2021.2} yields an ICE for passing complex components to `random_number()` or `log()`.<br>
    !>  This appears to have been fixed in \ifort{2021.4}.<br>
    !>  \remedy{1.5}
    !>  For now, the implementations are kept separately, since the installation of the new \ifort versions on supercomputers often lags.<br>
    !>  \remedy{Oct, 2022}
    !>  The `complex` interface of the routines is now deprecated and removed.<br>
    !>
    !>  \final{setExpRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setExpRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpRandDD_RK5(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandDD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(inout)                 :: rand
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpRandDD_RK4(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandDD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(inout)                 :: rand
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpRandDD_RK3(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandDD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(inout)                 :: rand
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpRandDD_RK2(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandDD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(inout)                 :: rand
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpRandDD_RK1(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandDD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(inout)                 :: rand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpRandSD_RK5(rand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpRandSD_RK4(rand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpRandSD_RK3(rand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpRandSD_RK2(rand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpRandSD_RK1(rand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setExpRandSM_RK5(rand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSM_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma, mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setExpRandSM_RK4(rand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSM_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma, mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setExpRandSM_RK3(rand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSM_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma, mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setExpRandSM_RK2(rand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSM_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma, mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setExpRandSM_RK1(rand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExpRandSM_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(inout)                 :: rand
        real(RKC)   , intent(in)                    :: sigma, mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distExp