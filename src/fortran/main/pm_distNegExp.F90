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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>Negative Exponential distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>Negative Exponential distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  The **PDF** of the **Negative Exponential distribution** with the two location and scale parameters \f$(\mu, \sigma)\f$ is defined as,<br>
!>  \f{equation}{
!>      \large
!>      F(x | \mu, \sigma) =
!>      \begin{cases}
!>          \frac{1}{\sigma} \exp(\frac{x - \mu}{\sigma})   &,~ -\infty < x \leq \mu \\
!>          0                                               &,~ x > \mu
!>      \end{cases}
!>  \f}
!>  where \f$-\infty < \mu < +\infty\f$ is the location parameter and \f$0 < \sigma < +\infty\f$ is the scale parameter of the distribution.<br>
!>
!>  The corresponding **CDF** with the two location and scale parameters \f$(\mu, \sigma)\f$ is defined as,<br>
!>  \f{equation}{
!>      \large
!>      \ms{CDF}(x | \mu, \sigma) =
!>      \begin{cases}
!>          \exp(\frac{x - \mu}{\sigma})    &,~ -\infty < x \leq \mu \\
!>          0                               &,~ x > \mu
!>      \end{cases}
!>  \f}
!>  where \f$-\infty < \mu < +\infty\f$ is the location parameter and \f$0 < \sigma < +\infty\f$ is the scale parameter of the distribution.<br>
!>
!>  **Random Number Generation**<br>
!>
!>  Assuming \f$U \in (0, 1]\f$ is a uniformly-distributed random variate,<br>
!>
!>  \f{equation}{
!>      \large
!>      T = \sigma\log(U) + \mu ~,
!>  \f}
!>
!>  is a random variate from the Negative Exponential distribution with the location parameter \f$-\infty < \mu < +\infty\f$
!>  and the scale parameter \f$0 < \sigma < +\infty\f$ such that \f$\mu \leq T < +\infty\f$.<br>
!>
!>  \note
!>  Note that the mean of the Negative Exponential distribution is the scale parameter: \f$\mathrm{mean} = \sigma\f$).<br>
!>
!>  \see
!>  [pm_distExp](@ref pm_distExp)<br>
!>  [pm_distUnif](@ref pm_distUnif)<br>
!>
!>  \test
!>  [test_pm_distNegExp](@ref test_pm_distNegExp)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distNegExp

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distNegExp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Negative Exponential
    !>  as defined in the description of [pm_distNegExp](@ref pm_distNegExp).
    !>
    !>  \details
    !>  See the documentation of [pm_distNegExp](@ref pm_distNegExp) for the definition of the Negative Exponential distribution.
    !>
    !>  \interface{distNegExp_type}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: distNegExp_type
    !>      type(distNegExp_type) :: distNegExp
    !>
    !>      distNegExp = distNegExp_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distNegExp_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distNegExp_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the Probability Density Function (PDF) of the
    !>  <b>Negative Exponential distribution</b> for an input `x` within the support of the distribution \f$(-\infty, \mu]\f$.
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
    !>  \interface{getNegExpLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: getNegExpLogPDF
    !>
    !>      logPDF = getNegExpLogPDF(x, mu = mu, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `x <= mu` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setNegExpLogPDF](@ref pm_distNegExp::setNegExpLogPDF)<br>
    !>
    !>  \example{getNegExpLogPDF}
    !>  \include{lineno} example/pm_distNegExp/getNegExpLogPDF/main.F90
    !>  \compilef{getNegExpLogPDF}
    !>  \output{getNegExpLogPDF}
    !>  \include{lineno} example/pm_distNegExp/getNegExpLogPDF/main.out.F90
    !>  \postproc{getNegExpLogPDF}
    !>  \include{lineno} example/pm_distNegExp/getNegExpLogPDF/main.py
    !>  \vis{getNegExpLogPDF}
    !>  \image html pm_distNegExp/getNegExpLogPDF/getNegExpLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getNegExpLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNegExpLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNegExpLogPDFXMI_RK5(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpLogPDFXMI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNegExpLogPDFXMI_RK4(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpLogPDFXMI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNegExpLogPDFXMI_RK3(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpLogPDFXMI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNegExpLogPDFXMI_RK2(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpLogPDFXMI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNegExpLogPDFXMI_RK1(x, mu, invSigma) result(logPDF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpLogPDFXMI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: logPDF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the Probability Density Function (PDF) of the
    !>  <b>Negative Exponential distribution</b> for an input `x` within the support of the distribution \f$(-\infty, \mu]\f$.
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
    !>  \interface{setNegExpLogPDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: setNegExpLogPDF
    !>
    !>      call setNegExpLogPDF(logPDF, x)
    !>      call setNegExpLogPDF(logPDF, x, mu)
    !>      call setNegExpLogPDF(logPDF, x, invSigma, logInvSigma)
    !>      call setNegExpLogPDF(logPDF, x, mu, invSigma, logInvSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `x <= mu` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  The condition `logInvSigma = log(invSigma)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getNegExpLogPDF](@ref pm_distNegExp::getNegExpLogPDF)<br>
    !>
    !>  \example{setNegExpLogPDF}
    !>  \include{lineno} example/pm_distNegExp/setNegExpLogPDF/main.F90
    !>  \compilef{setNegExpLogPDF}
    !>  \output{setNegExpLogPDF}
    !>  \include{lineno} example/pm_distNegExp/setNegExpLogPDF/main.out.F90
    !>  \postproc{setNegExpLogPDF}
    !>  \include{lineno} example/pm_distNegExp/setNegExpLogPDF/main.py
    !>  \vis{setNegExpLogPDF}
    !>  \image html pm_distNegExp/setNegExpLogPDF/setNegExpLogPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setNegExpLogPDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNegExpLogPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDDD_RK5(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDDD_RK4(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDDD_RK3(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDDD_RK2(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDDD_RK1(logPDF, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpLogPDFXMDD_RK5(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFXMDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpLogPDFXMDD_RK4(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFXMDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpLogPDFXMDD_RK3(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFXMDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpLogPDFXMDD_RK2(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFXMDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpLogPDFXMDD_RK1(logPDF, x, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFXMDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDIL_RK5(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDIL_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDIL_RK4(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDIL_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDIL_RK3(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDIL_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDIL_RK2(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDIL_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpLogPDFDIL_RK1(logPDF, x, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFDIL_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, invSigma, logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpLogPDFMIL_RK5(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFMIL_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpLogPDFMIL_RK4(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFMIL_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpLogPDFMIL_RK3(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFMIL_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpLogPDFMIL_RK2(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFMIL_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpLogPDFMIL_RK1(logPDF, x, mu, invSigma, logInvSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpLogPDFMIL_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: logPDF
        real(RKG)   , intent(in)                    :: x, mu, invSigma, logInvSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  <b>Negative Exponential distribution</b> for an input `x` within the support of the distribution \f$(-\infty, \mu]\f$.
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
    !>  \interface{getNegsetNegExpCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: getNegsetNegExpCDF
    !>
    !>      cdf = getNegExpCDF(x, mu = mu, invSigma = invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `x <= mu` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setNegExpCDF](@ref pm_distNegExp::setNegExpCDF)<br>
    !>
    !>  \example{getNegsetNegExpCDF}
    !>  \include{lineno} example/pm_distNegExp/getNegExpCDF/main.F90
    !>  \compilef{getNegsetNegExpCDF}
    !>  \output{getNegsetNegExpCDF}
    !>  \include{lineno} example/pm_distNegExp/getNegExpCDF/main.out.F90
    !>  \postproc{getNegsetNegExpCDF}
    !>  \include{lineno} example/pm_distNegExp/getNegExpCDF/main.py
    !>  \vis{getNegsetNegExpCDF}
    !>  \image html pm_distNegExp/getNegExpCDF/getNegExpCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{getNegsetNegExpCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNegExpCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getNegExpCDFXMI_RK5(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpCDFXMI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getNegExpCDFXMI_RK4(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpCDFXMI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getNegExpCDFXMI_RK3(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpCDFXMI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getNegExpCDFXMI_RK2(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpCDFXMI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getNegExpCDFXMI_RK1(x, mu, invSigma) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpCDFXMI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: x
        real(RKG)   , intent(in)    , optional  :: mu, invSigma
        real(RKG)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  <b>Negative Exponential distribution</b> for an input `x` within the support of the distribution \f$(-\infty, \mu]\f$.
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
    !>                              (**optional**, default = `0`. It can be present <b>if and only if</b> `invSigma` is also present.)
    !>  \param[in]  invSigma    :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `x`
    !>                              containing the inverse scale (i.e., the rate or the inverse mean) parameter of the distribution.<br>
    !>                              (**optional**, default = `1.`.)
    !>
    !>  \interface{setNegExpCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: setNegExpCDF
    !>
    !>      call setNegExpCDF(cdf, x)
    !>      call setNegExpCDF(cdf, x, invSigma)
    !>      call setNegExpCDF(cdf, x, mu, invSigma)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `x <= mu` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invSigma` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getNegExpCDF](@ref pm_distNegExp::getNegExpCDF)<br>
    !>
    !>  \example{setNegExpCDF}
    !>  \include{lineno} example/pm_distNegExp/setNegExpCDF/main.F90
    !>  \compilef{setNegExpCDF}
    !>  \output{setNegExpCDF}
    !>  \include{lineno} example/pm_distNegExp/setNegExpCDF/main.out.F90
    !>  \postproc{setNegExpCDF}
    !>  \include{lineno} example/pm_distNegExp/setNegExpCDF/main.py
    !>  \vis{setNegExpCDF}
    !>  \image html pm_distNegExp/setNegExpCDF/setNegExpCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \final{setNegExpCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNegExpCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpCDFXDD_RK5(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpCDFXDD_RK4(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpCDFXDD_RK3(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpCDFXDD_RK2(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpCDFXDD_RK1(cdf, x)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpCDFXDI_RK5(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpCDFXDI_RK4(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpCDFXDI_RK3(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpCDFXDI_RK2(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpCDFXDI_RK1(cdf, x, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXDI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpCDFXMI_RK5(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXMI_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpCDFXMI_RK4(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXMI_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpCDFXMI_RK3(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXMI_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpCDFXMI_RK2(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXMI_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpCDFXMI_RK1(cdf, x, mu, invSigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpCDFXMI_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: cdf
        real(RKG)   , intent(in)                    :: x, mu, invSigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Negative Exponential distribution,
    !>  optionally with the specified input location and scale parameters `mu, sigma`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distNegExp](@ref pm_distNegExp) for more details.<br>
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
    !>  \interface{getNegExpRand}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: getNegExpRand
    !>
    !>      rand = getNegExpRand(sigma, mu = mu)
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
    !>  [setNegExpRand](@ref pm_distNegExp::setNegExpRand)<br>
    !>
    !>  \example{getNegExpRand}
    !>  \include{lineno} example/pm_distNegExp/getNegExpRand/main.F90
    !>  \compilef{getNegExpRand}
    !>  \output{getNegExpRand}
    !>  \include{lineno} example/pm_distNegExp/getNegExpRand/main.out.F90
    !>  \postproc{getNegExpRand}
    !>  \include{lineno} example/pm_distNegExp/getNegExpRand/main.py
    !>  \vis{getNegExpRand}
    !>  \image html pm_distNegExp/getNegExpRand/getNegExpRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
    !>
    !>  \final{getNegExpRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getNegExpRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getNegExpRandSM_RK5(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpRandSM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                                   :: rand
        real(RKG)   , intent(in)                    :: sigma
        real(RKG)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getNegExpRandSM_RK4(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpRandSM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                                   :: rand
        real(RKG)   , intent(in)                    :: sigma
        real(RKG)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getNegExpRandSM_RK3(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpRandSM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                                   :: rand
        real(RKG)   , intent(in)                    :: sigma
        real(RKG)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getNegExpRandSM_RK2(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpRandSM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                                   :: rand
        real(RKG)   , intent(in)                    :: sigma
        real(RKG)   , intent(in)    , optional      :: mu
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getNegExpRandSM_RK1(sigma, mu) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNegExpRandSM_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                                   :: rand
        real(RKG)   , intent(in)                    :: sigma
        real(RKG)   , intent(in)    , optional      :: mu
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank of) random value(s) from the Negative Exponential distribution,
    !>  optionally with the specified input location and scale parameters `mu, sigma`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distNegExp](@ref pm_distNegExp) for more details.<br>
    !>
    !>  \param[out] rand    :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `urand`,
    !>                          containing random value(s) from the distribution.<br>
    !>  \param[in]  urand   :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                          <ul>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing uniformly-distributed random value(s) with the range `0 <= urand < 1`.<br>
    !>                          Such random value(s) can be readily obtained via the Fortran intrinsic procedure `random_number()` or via [getUnifRand()](@ref pm_distUnif::getUnifRand).<br>
    !>                          Supplying this argument ensures the purity of the procedures, allowing further compiler optimizations.<br>
    !>  \param[in]  sigma   :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `urand`,
    !>                          representing the scale parameter of the distribution.<br>
    !>                          (**optional**, default = `1.`.)
    !>  \param[in]  mu      :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of the same type and kind as `urand`,
    !>                          representing the location parameter of the distribution.<br>
    !>                          (**optional**, default = `0`. It can be present <b>if and only if</b> `sigma` is also present.)
    !>
    !>  \interface{setNegExpRand}
    !>  \code{.F90}
    !>
    !>      use pm_distNegExp, only: setNegExpRand
    !>
    !>      call setNegExpRand(rand, urand)
    !>      call setNegExpRand(rand, urand, sigma)
    !>      call setNegExpRand(rand, urand, sigma, mu)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < sigma` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= urand < 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getNegExpRand](@ref pm_distNegExp::getNegExpRand)<br>
    !>
    !>  \example{setNegExpRand}
    !>  \include{lineno} example/pm_distNegExp/setNegExpRand/main.F90
    !>  \compilef{setNegExpRand}
    !>  \output{setNegExpRand}
    !>  \include{lineno} example/pm_distNegExp/setNegExpRand/main.out.F90
    !>  \postproc{setNegExpRand}
    !>  \include{lineno} example/pm_distNegExp/setNegExpRand/main.py
    !>  \vis{setNegExpRand}
    !>  \image html pm_distNegExp/setNegExpRand/setNegExpRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distNegExp](@ref test_pm_distNegExp)
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
    !>  \final{setNegExpRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setNegExpRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpRandUDD_RK5(rand, urand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpRandUDD_RK4(rand, urand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpRandUDD_RK3(rand, urand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpRandUDD_RK2(rand, urand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpRandUDD_RK1(rand, urand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpRandUSD_RK5(rand, urand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpRandUSD_RK4(rand, urand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpRandUSD_RK3(rand, urand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpRandUSD_RK2(rand, urand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpRandUSD_RK1(rand, urand, sigma)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setNegExpRandUSM_RK5(rand, urand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma, mu
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setNegExpRandUSM_RK4(rand, urand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma, mu
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setNegExpRandUSM_RK3(rand, urand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma, mu
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setNegExpRandUSM_RK2(rand, urand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma, mu
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setNegExpRandUSM_RK1(rand, urand, sigma, mu)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNegExpRandUSM_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: rand
        real(RKG)   , intent(in)                    :: urand, sigma, mu
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distNegExp