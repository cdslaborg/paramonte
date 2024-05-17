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
!>  This module contains classes and procedures for computing various statistical quantities related to the <b>LogUniform (or Reciprocal) distribution</b>.
!>
!>  \details
!>  Specifically, this module contains routines for computing the following quantities of the <b>LogUniform (or Reciprocal) distribution</b>:<br>
!>  <ol>
!>      <li>    the Probability Density Function (**PDF**)
!>      <li>    the Cumulative Distribution Function (**CDF**)
!>      <li>    the Random Number Generation from the distribution (**RNG**)
!>      <li>    the Inverse Cumulative Distribution Function (**ICDF**) or the **Quantile Function**
!>  </ol>
!>
!>  In probability and statistics, the **reciprocal distribution**, also known as the **LogUniform distribution**, is a continuous probability distribution.<br>
!>  It is characterized by its probability density function, within the support of the distribution, being proportional to the reciprocal of the variable.<br>
!>
!>  The **PDF** of the <b>LogUniform distribution</b> with the two scale parameters \f$(x_\mathrm{min}, x_\mathrm{max})\f$
!>  over a strictly-positive support \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$ is defined as,<br>
!>
!>  \f{eqnarray}{
!>      \large
!>      \pi(x | x_\mathrm{min}, x_\mathrm{max})
!>      &=& \frac{\eta(x_\mathrm{min}, x_\mathrm{max})}{x} ~, \\
!>      &=& \frac{1}{x ~ \left[\log(x_\mathrm{min}) - \log(x_\mathrm{max})\right]} ~,
!>  \f}
!>
!>  where \f$\log(\cdot)\f$ is the natural logarithm and the condition \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ holds.<br>
!>  The term \f$\eta(x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distLogUnif::getLogUnifPDFNF) of the PDF,<br>
!>
!>  \f{equation}{
!>      \large
!>      \eta(x_\mathrm{min}, x_\mathrm{max}) = \frac{1}{\log(x_\mathrm{min}) - \log(x_\mathrm{max})} ~,
!>  \f}
!>
!>  The LogUniform distribution can be considered an special case of the [Truncated Pareto distribution](@ref pm_distPareto),
!>  or the [Truncated Power distribution](@ref pm_distPower), or the [Truncated Poweto distribution](@ref pm_distPoweto).<br>
!>
!>  The **CDF** of the LogUniform distribution is detailed in [pm_distLogUnif](@ref pm_distLogUnif).<br>
!>  The corresponding **CDF** of the LogUniform distribution is given by,
!>  \f{eqnarray}{
!>      \large
!>      \mathrm{CDF}(x | x_\mathrm{min}, x_\mathrm{max})
!>      &=& \eta(x_\mathrm{min}, x_\mathrm{max}) \log\left(\frac{x}{x_\mathrm{xmin}}\right) ~, \\
!>      &=& 1 + \eta(x_\mathrm{min}, x_\mathrm{max}) \log\left(\frac{x}{x_\mathrm{xmax}}\right) ~, \\
!>  \f}
!>  where \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ holds.<br>
!>  The term \f$\eta(x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distLogUnif::getLogUnifPDFNF) of the PDF,<br>
!>  \f{eqnarray}{
!>      \large
!>      \eta(x_\mathrm{min}, x_\mathrm{max})
!>      &=& \frac{1}{\log(x_\mathrm{max}) - \log(x_\mathrm{min})} ~.
!>  \f}
!>
!>  The corresponding **Inverse CDF** or **Quantile Function** of the LogUniform distribution is given by,
!>  \f{equation}{
!>      \large
!>      Q(\mathrm{CDF}(x); x_\mathrm{min}, x_\mathrm{max}) \equiv x = x_\mathrm{min} ~ \exp\left(\frac{\mathrm{CDF}(x)}{\eta(x_\mathrm{min}, x_\mathrm{max})}\right) ~,
!>  \f}
!>  where the condition \f$\mathbf{0 < x_\mathrm{min} \leq Q(\mathrm{CDF}(x); x_\mathrm{min}, x_\mathrm{max}) \leq x_\mathrm{max} < +\infty}\f$ holds
!>  and \f$\eta(x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distLogUnif::getLogUnifPDFNF) of the PDF,<br>
!>  \f{equation}{
!>      \large
!>      \eta(x_\mathrm{min}, x_\mathrm{max}) = \frac{1}{\log(x_\mathrm{min}) - \log(x_\mathrm{max})} ~,
!>  \f}
!>
!>  <b>Random Number Generation</b><br>
!>
!>  Assuming that \f$U \in [0, 1]\f$ is a uniformly-distributed random variate, the transformed random variable,
!>  \f{equation}{
!>      \large
!>      x = x_\mathrm{min} ~ \exp\left(\frac{U}{\eta(x_\mathrm{min}, x_\mathrm{max})}\right) ~,
!>  \f}
!>  follows a **LogUniform distribution** with parameters \f$(x_\mathrm{min}, x_\mathrm{max}\f$)
!>  where \f$\mathbf{0 < x_\mathrm{min} \leq x \leq x_\mathrm{max} < +\infty}\f$ holds
!>  and \f$\eta(x_\mathrm{min}, x_\mathrm{max})\f$ is the [normalization factor](@ref pm_distLogUnif::getLogUnifPDFNF) of the PDF,<br>
!>  \f{equation}{
!>      \large
!>      \eta(x_\mathrm{min}, x_\mathrm{max}) = \frac{1}{\log(x_\mathrm{min}) - \log(x_\mathrm{max})} ~,
!>  \f}
!>
!>  \see
!>  [pm_distLogUnif](@ref pm_distLogUnif)<br>
!>  [pm_distLogUnif](@ref pm_distLogUnif)<br>
!>  [pm_distLogUnif](@ref pm_distLogUnif)<br>
!>  [pm_distPareto](@ref pm_distPareto)<br>
!>  [pm_distPoweto](@ref pm_distPoweto)<br>
!>  [pm_distPower](@ref pm_distPower)<br>
!>
!>  \test
!>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distLogUnif

    use pm_kind, only: SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distLogUnif"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type LogUniform
    !>  as defined in the description of [pm_distLogUnif](@ref pm_distLogUnif).
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for the definition of the LogUniform distribution.
    !>
    !>  \interface{distLogUnif_type}
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: distLogUnif_type
    !>      type(distLogUnif_type) :: distLogUnif
    !>
    !>      distLogUnif = distLogUnif_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \finmain{distLogUnif_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distLogUnif_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the normalization factor of the PDF
    !>  of the LogUniform distribution for an input parameter set \f$(x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \brief
    !>  The normalization factor \f$\eta(x_\mathrm{min}, x_\mathrm{max})\f$ of the PDF of the LogUniform distribution is,
    !>
    !>  \f{equation}{
    !>
    !>      \large
    !>      \eta(x_\mathrm{min}, x_\mathrm{max}) = \frac{1}{\log(x_\mathrm{max}) - \log(x_\mathrm{min})} ~,~ 0 < x_\mathrm{min} < x_\mathrm{max} < +\infty
    !>
    !>  \f}
    !>
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.<br>
    !>
    !>  The primary use of this interface is to compute the normalization factor of the PDF of the LogUniform distribution for a fixed set of parameters
    !>  and use it in subsequent repeated calculations of the properties of the LogUniform distribution to improve the runtime performance by eliminating redundant calculations.<br>
    !>
    !>  \param[in]  logMinX :   The input scalar, or array of the same shape as other array-like arguments, of type `real` of kind \RKALL,
    !>                          containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logMinX`,
    !>                          containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `pdfnf`             :   The output scalar or array of the same shape as any input array-like argument,
    !>                          of the same type and kind as the input argument `logMinX`, containing
    !>                          the normalization factor of the PDF of LogUniform distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: getLogUnifPDFNF
    !>
    !>      pdfnf = getLogUnifPDFNF(logMinX, logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogUnifPDF](@ref pm_distLogUnif::getLogUnifPDF)<br>
    !>  [setLogUnifPDF](@ref pm_distLogUnif::setLogUnifPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifPDFNF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifPDFNF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifPDFNF/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/getLogUnifPDFNF/getLogUnifPDFNF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogUnifPDFNF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogUnifPDFNF_RK5(logMinX, logMaxX) result(pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFNF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logMinX, logMaxX
        real(RKC)                               :: pdfnf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogUnifPDFNF_RK4(logMinX, logMaxX) result(pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFNF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logMinX, logMaxX
        real(RKC)                               :: pdfnf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogUnifPDFNF_RK3(logMinX, logMaxX) result(pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFNF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logMinX, logMaxX
        real(RKC)                               :: pdfnf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogUnifPDFNF_RK2(logMinX, logMaxX) result(pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFNF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logMinX, logMaxX
        real(RKC)                               :: pdfnf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogUnifPDFNF_RK1(logMinX, logMaxX) result(pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFNF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logMinX, logMaxX
        real(RKC)                               :: pdfnf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Probability Density Function (PDF) of the
    !>  LogUniform distribution for an input `x` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.<br>
    !>
    !>  \param[in]  x       :   The input scalar or array of the same shape as other array like arguments,
    !>                          of type `real` of kind \RKALL, containing the value at which the PDF must be computed.<br>
    !>  \param[in]  minx    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                          containing the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  maxx    :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                          containing the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `pdf`               :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                          and kind the input argument `x`, containing the PDF of the distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: getLogUnifPDF
    !>
    !>      pdf = getLogUnifPDF(x, minx, maxx)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. < minx` must hold for the corresponding input arguments.<br>
    !>  The condition `minx <= x` must hold for the corresponding input arguments.<br>
    !>  The condition `x <= maxx` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogUnifPDF](@ref pm_distLogUnif::setLogUnifPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifPDF/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/getLogUnifPDF/getLogUnifPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogUnifPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogUnifPDFMM_RK5(x, minx, maxx) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFMM_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: x, minx, maxx
        real(RKC)                               :: pdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogUnifPDFMM_RK4(x, minx, maxx) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFMM_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: x, minx, maxx
        real(RKC)                               :: pdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogUnifPDFMM_RK3(x, minx, maxx) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFMM_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: x, minx, maxx
        real(RKC)                               :: pdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogUnifPDFMM_RK2(x, minx, maxx) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFMM_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: x, minx, maxx
        real(RKC)                               :: pdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogUnifPDFMM_RK1(x, minx, maxx) result(pdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifPDFMM_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: x, minx, maxx
        real(RKC)                               :: pdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Probability Density Function (PDF) of the
    !>  LogUniform distribution for an input `x` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.<br>
    !>
    !>  \param[out] pdf     :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                          and kind the input argument `x`, containing the natural logarithm of the PDF of the distribution.
    !>  \param[in]  x       :   The input scalar or array of the same shape as other array like arguments,
    !>                          of type `real` of kind \RKALL, containing the value at which the PDF must be computed.<br>
    !>  \param[in]  pdfnf   :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `x`,
    !>                          containing the normalization factor of the PDF of the LogUniform distribution.<br>
    !>                          Specifying this argument when calling this procedure repeatedly with fixed \f$(x_\mathrm{min}, x_\mathrm{max})\f$
    !>                          parameters will significantly improve the runtime performance.<br>
    !>                          This argument can be readily obtained by calling [getLogUnifPDFNF(logMinX, logMaxX)](@ref pm_distLogUnif::getLogUnifPDFNF).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: setLogUnifPDF
    !>
    !>      call setLogUnifPDF(pdf, x, pdfnf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Be warned that the procedures under this generic interface have no mechanism for checking for the consistency of the specified input arguments with each other.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogUnifPDF](@ref pm_distLogUnif::setLogUnifPDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifPDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifPDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifPDF/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/setLogUnifPDF/setLogUnifPDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setLogUnifPDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogUnifPDF_RK5(pdf, x, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifPDF_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, pdfnf
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogUnifPDF_RK4(pdf, x, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifPDF_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, pdfnf
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogUnifPDF_RK3(pdf, x, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifPDF_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, pdfnf
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogUnifPDF_RK2(pdf, x, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifPDF_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, pdfnf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogUnifPDF_RK1(pdf, x, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifPDF_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: pdf
        real(RKC)   , intent(in)                    :: x, pdfnf
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Cumulative Distribution Function (CDF) of the
    !>  LogUniform distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.
    !>
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the natural logarithm of the value at which the CDF must be computed.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `cdf`                   :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `logx`, containing the CDF of the distribution.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: getLogUnifCDF
    !>
    !>      cdf = getLogUnifCDF(logx, logMinX, logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  The conditions `logx <= logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogUnifCDF](@ref pm_distLogUnif::setLogUnifCDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifCDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifCDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifCDF/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/getLogUnifCDF/getLogUnifCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogUnifCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogUnifCDFLL_RK5(logx, logMinX, logMaxX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifCDFLL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: logx, logMinX, logMaxX
        real(RKC)                               :: cdf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogUnifCDFLL_RK4(logx, logMinX, logMaxX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifCDFLL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: logx, logMinX, logMaxX
        real(RKC)                               :: cdf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogUnifCDFLL_RK3(logx, logMinX, logMaxX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifCDFLL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: logx, logMinX, logMaxX
        real(RKC)                               :: cdf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogUnifCDFLL_RK2(logx, logMinX, logMaxX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifCDFLL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: logx, logMinX, logMaxX
        real(RKC)                               :: cdf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogUnifCDFLL_RK1(logx, logMinX, logMaxX) result(cdf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifCDFLL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: logx, logMinX, logMaxX
        real(RKC)                               :: cdf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Cumulative Distribution Function (CDF) of the
    !>  LogUniform distribution for an input `log(x)` within the support of the distribution \f$x \in [0 < x_\mathrm{min}, x_\mathrm{max} < +\infty]\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.
    !>
    !>  \param[out] cdf         :   The output scalar or array of the same shape as any input array-like argument, of the same type
    !>                              and kind the input argument `x`, containing the CDF of the distribution.
    !>  \param[in]  logx        :   The input scalar or array of the same shape as other array like arguments, of type `real` of kind \RKALL,
    !>                              containing the natural logarithm of the value at which the CDF must be computed.<br>
    !>  \param[in]  logMinX     :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  pdfnf       :   The input scalar or array of the same shape as other array-like arguments, of the same type and kind as `logx`,
    !>                              containing the normalization factor of the PDF of the LogUniform distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getLogUnifPDFNF(logMinX, logMaxX)](@ref pm_distLogUnif::getLogUnifPDFNF).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: setLogUnifCDF
    !>
    !>      call setLogUnifCDF(cdf, logx, logMinX, pdfnf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `logMinX <= logx` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Be warned that the procedures under this generic interface have no mechanism for checking for the consistency of the specified input arguments with each other.<br>
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogUnifCDF](@ref pm_distLogUnif::setLogUnifCDF)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifCDF/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifCDF/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifCDF/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/setLogUnifCDF/setLogUnifCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setLogUnifCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogUnifCDFLL_RK5(cdf, logx, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifCDFLL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: logx, logMinX, pdfnf
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogUnifCDFLL_RK4(cdf, logx, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifCDFLL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: logx, logMinX, pdfnf
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogUnifCDFLL_RK3(cdf, logx, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifCDFLL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: logx, logMinX, pdfnf
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogUnifCDFLL_RK2(cdf, logx, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifCDFLL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: logx, logMinX, pdfnf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogUnifCDFLL_RK1(cdf, logx, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifCDFLL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: cdf
        real(RKC)   , intent(in)                    :: logx, logMinX, pdfnf
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of the natural logarithm(s) of quantile corresponding to
    !>  the specified CDF of <b>LogUniform distribution</b> with parameters \f$(x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the quantile of the LogUniform distribution.<br>
    !>
    !>  \param[in]  cdf         :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `logMinX`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  logMaxX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `logMinX`,
    !>                              containing the natural logarithm of the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `logx`                  :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `logMinX`, containing the natural logarithm of the quantile corresponding to the input `cdf`.<br>
    !>                              By definition, the condition `logMinX <= logx <= logMaxX` holds.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: getLogUnifLogQuan
    !>
    !>      logx = getLogUnifLogQuan(cdf, logMinX, logMaxX)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= cdf <= 1` must hold for the corresponding input arguments.<br>
    !>  The conditions `logMinX < logMaxX` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogUnifLogQuan](@ref pm_distLogUnif::setLogUnifLogQuan)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifLogQuan/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifLogQuan/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifLogQuan/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/getLogUnifLogQuan/getLogUnifLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogUnifLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogUnifLogQuanLL_RK5(cdf, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifLogQuanLL_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: cdf, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogUnifLogQuanLL_RK4(cdf, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifLogQuanLL_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: cdf, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogUnifLogQuanLL_RK3(cdf, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifLogQuanLL_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: cdf, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogUnifLogQuanLL_RK2(cdf, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifLogQuanLL_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: cdf, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogUnifLogQuanLL_RK1(cdf, logMinX, logMaxX) result(logx)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifLogQuanLL_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: cdf, logMinX, logMaxX
        real(RKC)                               :: logx
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the natural logarithm(s) of quantile corresponding to
    !>  the specified CDF of <b>LogUniform distribution</b> with parameters \f$(x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the quantile of the LogUniform distribution.<br>
    !>
    !>  \param[out] logx        :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                              of the same type and kind as `cdf`, containing the natural logarithm of the quantile corresponding to the input `cdf`.<br>
    !>                              By definition, the condition `logMinX <= logx <= logMaxX` holds.
    !>  \param[in]  cdf         :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of
    !>                              <ul>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the natural logarithm of the desired CDF value of the distribution corresponding to the output quantile.<br>
    !>  \param[in]  logMinX     :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `cdf`,
    !>                              containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  pdfnf       :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `cdf`,
    !>                              containing the natural logarithm of the normalization factor of the PDF of the LogUniform distribution.<br>
    !>                              Specifying this argument when calling this procedure repeatedly with fixed \f$(x_\mathrm{min}, x_\mathrm{max})\f$
    !>                              parameters will significantly improve the runtime performance.<br>
    !>                              This argument can be readily obtained by calling [getLogUnifPDFNF(logMinX, logMaxX)](@ref pm_distLogUnif::getLogUnifPDFNF).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: setLogUnifLogQuan
    !>
    !>      call setLogUnifLogQuan(logx, cdf, logMinX, pdfnf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= cdf <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogUnifLogQuan](@ref pm_distLogUnif::getLogUnifLogQuan)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifLogQuan/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifLogQuan/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifLogQuan/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/setLogUnifLogQuan/setLogUnifLogQuan.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setLogUnifLogQuan

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogUnifLogQuanLLLP_RK5(logx, cdf, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogQuanLLLP_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: cdf, logMinX, pdfnf
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogUnifLogQuanLLLP_RK4(logx, cdf, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogQuanLLLP_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: cdf, logMinX, pdfnf
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogUnifLogQuanLLLP_RK3(logx, cdf, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogQuanLLLP_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: cdf, logMinX, pdfnf
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogUnifLogQuanLLLP_RK2(logx, cdf, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogQuanLLLP_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: cdf, logMinX, pdfnf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogUnifLogQuanLLLP_RK1(logx, cdf, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogQuanLLLP_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logx
        real(RKC)   , intent(in)                    :: cdf, logMinX, pdfnf
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or array of arbitrary rank) of random value(s) from the
    !>  <b>LogUniform distribution</b> with parameters \f$(x_\mathrm{min}, x_\mathrm{max})\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.<br>
    !>
    !>  \param[in]  minx    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of
    !>                          <ul>
    !>                              <li>    type `integer` of kind \IKALL, <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  maxx    :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `minx`,
    !>                          containing the second scale parameter of the distribution, representing the maximum of the support of the distribution.<br>
    !>
    !>  \return
    !>  `rand`              :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                          of the same type and kind as `minx`, containing the random value(s) from the specified distribution.<br>
    !>                          By definition, the condition `minx <= rand < maxx` holds.<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: getLogUnifRand
    !>
    !>      rand = getLogUnifRand(minx, maxx)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `minx < maxx` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setLogUnifLogRand](@ref pm_distLogUnif::setLogUnifLogRand)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifRand/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/getLogUnifRand/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/getLogUnifRand/getLogUnifRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to `complex` arguments.<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogUnifRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    impure elemental module function getLogUnifRandMM_IK5(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC), intent(in)                :: minx, maxx
        integer(IKC)                            :: rand
    end function
#endif

#if IK4_ENABLED
    impure elemental module function getLogUnifRandMM_IK4(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC), intent(in)                :: minx, maxx
        integer(IKC)                            :: rand
    end function
#endif

#if IK3_ENABLED
    impure elemental module function getLogUnifRandMM_IK3(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC), intent(in)                :: minx, maxx
        integer(IKC)                            :: rand
    end function
#endif

#if IK2_ENABLED
    impure elemental module function getLogUnifRandMM_IK2(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC), intent(in)                :: minx, maxx
        integer(IKC)                            :: rand
    end function
#endif

#if IK1_ENABLED
    impure elemental module function getLogUnifRandMM_IK1(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC), intent(in)                :: minx, maxx
        integer(IKC)                            :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getLogUnifRandMM_RK5(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                :: minx, maxx
        real(RKC)                               :: rand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getLogUnifRandMM_RK4(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                :: minx, maxx
        real(RKC)                               :: rand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getLogUnifRandMM_RK3(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                :: minx, maxx
        real(RKC)                               :: rand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getLogUnifRandMM_RK2(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                :: minx, maxx
        real(RKC)                               :: rand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getLogUnifRandMM_RK1(minx, maxx) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogUnifRandMM_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                :: minx, maxx
        real(RKC)                               :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar (or array of arbitrary rank) of the natural logarithm(s) of random value(s) from the
    !>  <b>LogUniform distribution</b> with parameters \f$(x_\mathrm{min}, x_\mathrm{max})\f$.
    !>
    !>  \details
    !>  See the documentation of [pm_distLogUnif](@ref pm_distLogUnif) for more information on the LogUniform distribution.<br>
    !>
    !>  \param[out] logRand :   The output scalar (or array of the same rank, shape, and size as other array like arguments),
    !>                          of the same type and kind as `urand`, containing the natural logarithm of the randomly-generated value(s) from the specified distribution.<br>
    !>                          By definition, the condition `logMinX <= logRand < logMaxX` holds.
    !>  \param[in]  urand   :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of
    !>                          <ul>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing a random value from the standard Uniform distribution with support `[0, 1)`.<br>
    !>                          This argument can be readily obtained by calling [getUnifRand(0., 1.)](@ref pm_distUnif::getUnifRand) for the desired `real` kind.<br>
    !>  \param[in]  logMinX :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `urand`,
    !>                          containing the natural logarithm of the first scale parameter of the distribution, representing the minimum of the support of the distribution.<br>
    !>  \param[in]  pdfnf   :   The input scalar (or array of the same rank, shape, and size as other array like arguments), of the same type and kind as `urand`,
    !>                          containing the natural logarithm of the normalization factor of the PDF of the LogUniform distribution.<br>
    !>                          Specifying this argument when calling this procedure repeatedly with fixed \f$(x_\mathrm{min}, x_\mathrm{max})\f$
    !>                          parameters will significantly improve the runtime performance.<br>
    !>                          This argument can be readily obtained by calling [getLogUnifPDFNF(logMinX, logMaxX)](@ref pm_distLogUnif::getLogUnifPDFNF).<br>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_distLogUnif, only: setLogUnifLogRand
    !>
    !>      call setLogUnifLogRand(logRand, urand, logMinX, pdfnf)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= urand <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogUnifRand](@ref pm_distLogUnif::getLogUnifRand)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifLogRand/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifLogRand/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_distLogUnif/setLogUnifLogRand/main.py
    !>  \vis
    !>  \image html pm_distLogUnif/setLogUnifLogRand/setLogUnifLogRand.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distLogUnif](@ref test_pm_distLogUnif)
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setLogUnifLogRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setLogUnifLogRandLLLP_RK5(logRand, urand, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogRandLLLP_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: urand, logMinX, pdfnf
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setLogUnifLogRandLLLP_RK4(logRand, urand, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogRandLLLP_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: urand, logMinX, pdfnf
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setLogUnifLogRandLLLP_RK3(logRand, urand, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogRandLLLP_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: urand, logMinX, pdfnf
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setLogUnifLogRandLLLP_RK2(logRand, urand, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogRandLLLP_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: urand, logMinX, pdfnf
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setLogUnifLogRandLLLP_RK1(logRand, urand, logMinX, pdfnf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogUnifLogRandLLLP_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: logRand
        real(RKC)   , intent(in)                    :: urand, logMinX, pdfnf
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distLogUnif