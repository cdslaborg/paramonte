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
!>  This module contains procedures and generic interfaces for computing the Band photon distribution
!>  widely used in modeling the spectra of a class of celestial objects known as [Gamma-Ray Bursts](https://en.wikipedia.org/wiki/Gamma-ray_burst).<br>
!>
!>  \details
!>  The Band model is an empirical spectral model most widely used to fit GRB spectra,
!>  first proposed in *Band et al. 1993,  BATSE Observations of Gamma-Ray Burst Spectra. I. Spectral Diversity*.<br>
!>
!>  The model is **continuously differentiable** (i.e., its derivative is a continuous function).<br>
!>  It is characterized by four parameters:<br>
!>  <ol>
!>      <li>    the amplitude \f$A\f$ in \f$\ms{photons} \ms{s}^{-1} \ms{cm}^{-2} \kev^{-1}\f$,<br>
!>      <li>    a low-energy spectral index \f$\alpha\f$,<br>
!>      <li>    a high-energy spectral index \f$\beta < \alpha\f$,<br>
!>      <li>    a \f$\nu F_\nu\f$ peak energy \f$\epeak\f$ in units of \f$\kev\f$ (or equivalently, the **break energy**, \f$\ebreak\f$).<br>
!>  </ol>
!>  The \f$\nu F_\nu\f$ is the **photon spectrum** \f$f(E)\f$ integrated twice over all energies (\f$E^2f(E)\f$).<br>
!>  Therefore, \f$\nu F_\nu\f$ represents the total energy flux per energy band (i.e., **power density spectrum**).<br>
!>  The \f$\alpha\f$ index characterizes an asymptotic power-law (i.e., the tangential slope determined at \f$E\rightarrow 0\f$ in a logarithmic scale).<br>
!>  This may not characterize the actual low-energy power-law, determined within the data energy range when the **e-folding energy** denoted by \f$E_{0}\f$ approaches the lower energy bound.<br>
!>  Although the model was originally constructed based on the observed time-integrated BATSE catalog spectra, it has now become common practice to use the model to fit time-resolved GRB spectra as well.<br>
!>  There are, however, some time-resolved spectra that cannot be adequately fitted with this model.<br>
!>  The Band model has the following mathematical form,
!>  \f{equation}{
!>      \large
!>      f_{\ms{BAND}}(E) =
!>      \begin{cases}
!>          A\left(\frac{E}{100\kev}\right)^\alpha \exp\left(-\frac{\alpha - \beta}{\ebreak}E\right) &,~ \ms{if} & E < \ebreak \\
!>          A\left[\frac{\ebreak}{100\kev}\right]^{\alpha - \beta} \exp\left(\beta - \alpha\right) \left(\frac{E}{100\kev}\right)^\beta &,~ \ms{if} & E \geq \ebreak \\
!>      \end{cases}
!>  \f}
!>  where \f$\ebreak = \frac{\alpha - \beta}{2 + \alpha}\epeak = (\alpha - \beta)\efold\f$ is called the **break energy** of the model.<br>
!>
!>  The above formulation takes a **unit-full** input value for the energy, \f$E\f$, at which the spectrum must be computed.<br>
!>  Assuming the normalization energy is unity (i.e., \f$1\kev\f$ instead of \f$100\kev\f$ without loss of generality),
!>  the **normalized bounded unit-less** formulation, corresponding to the <b>probability density function (PDF)</b> for arbitrary \f$(\alpha, \beta)\f$takes the form,<br>
!>  \f{equation}{
!>      \large
!>      \pi_{\ms{BAND}}(x | \alpha, \beta, \xbreak) = \eta(\alpha, \beta, \xbreak)
!>      \begin{cases}
!>          x^\alpha \exp\left(-\frac{x}{\sigma}\right) &,~ \ms{if} & 0 < \ms{lb} \leq x < \xbreak \\
!>          \zeta x^\beta &,~ \ms{if} & \xbreak \leq x < \ms{ub} < +\infty \\
!>      \end{cases}
!>  \f}
!>  where,
!>  <ol>
!>      <li>    \f$x = \frac{E}{\kev}\f$,
!>      <li>    \f$\xbreak = \frac{\ebreak}{\kev}\f$,
!>      <li>    \f$\sigma = \frac{\efold}{\kev} = \frac{(\alpha - \beta)\kev}{\xbreak}\f$ is the **normalized e-folding energy**,
!>              a measure of the **scale** of the distribution below and above which the distribution
!>              approaches power-law behavior with exponents \f$\alpha\f$ and \f$\beta\f$ respectively,
!>      <li>    the factor \f$\eta(\alpha, \beta, \xbreak)\f$ is a normalization constant that properly normalizes the PDF,
!>      <li>    the factor \f$\zeta = \xbreak^{\alpha - \beta} \exp\left(\beta - \alpha\right)\f$ is
!>              the **coefficient of continuity** of the distribution that makes the distribution **continuously differentiable**.<br>
!>      <li>    the constants \f$(\ms{lb}, \ms{ub})\f$ represent the lower and upper bounds of the PDF, respectively.<br>
!>  </ol>
!>  When the condition \f$\alpha > -1\f$ holds, the lower component of the distribution
!>  follows the mathematical form of the PDF of the [Gamma distribution](@ref pm_distGamma),
!>  \f{eqnarray}{
!>      \large
!>      x^\alpha \exp\left(-\frac{x}{\sigma}\right)
!>      &=& \sigma^\alpha \left(\frac{x}{\sigma}\right)^\alpha \exp\left(-\frac{x}{\sigma}\right) \\
!>      &=& \sigma^{\alpha + 1}\Gamma(\alpha + 1)\pi_\mathcal{G}(x | \kappa = \alpha + 1, \sigma) ~,
!>  \f}
!>
!>  \see
!>  Kaneko, 2005, Spectral studies of gamma-ray burst prompt emission.<br>
!>  Band et al. 1993,  BATSE Observations of Gamma-Ray Burst Spectra. I. Spectral Diversity.<br>
!>  Eqn. A6 in [Shahmoradi and Nemiroff, 2015, MNRAS, 451:4645-4662.](https://www.cdslab.org/pubs/2015_Shahmoradi_I.pdf)<br>
!>
!>  \test
!>  [test_pm_distBand](@ref test_pm_distBand)<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distBand

    use pm_kind, only: SK, IK, LK, RKB

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distBand"

    !>  \brief
    !>  The scalar constant of type `real` of kind \RKB,
    !>  containing the average reported value for the low-energy exponent of the Band photon distribution model \f$\alpha\f$.<br>
    !>  This reported average value is primarily based on the analyses of data collected by
    !>  the BATSE telescope onboard the Compton Gamma-Ray Observatory (CGRO).<br>
    !>
    !>  \warning
    !>  This average value makes the Band distribution unnormalized because the integral of the PDF does not converge for \f$x\rightarrow0\f$.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    real(RKB), parameter :: MEAN_ALPHA = -1.1_RKB

    !>  \brief
    !>  The scalar constant of type `real` of kind \RKB,
    !>  containing the average reported value for the high-energy exponent of the Band photon distribution model \f$\beta\f$.<br>
    !>  This reported average value is primarily based on the analyses of data collected by
    !>  the BATSE telescope onboard the Compton Gamma-Ray Observatory (CGRO).<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA
    real(RKB), parameter :: MEAN_BETA = -2.3_RKB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for signifying distributions that are of type Band
    !>  as defined in the description of [pm_distBand](@ref pm_distBand).
    !>
    !>  \details
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for the definition of the Band distribution.
    !>
    !>  \interface{distBand_type}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: distBand_type
    !>      type(distBand_type) :: distBand
    !>
    !>      distBand = distBand_type()
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  This derived type is currently devoid of any components or type-bound procedures because of
    !>  the lack of portable and reliable support for Parameterized Derived Types (PDT) in some Fortran compilers.<br>
    !>  For now, the utility of this derived type is limited to generic interface resolutions.<br>
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \todo
    !>  \pvhigh
    !>  This derived type must be converted to PDT and the relevant components and methods must be added once PDTs are well supported.
    !>
    !>  \final{distBand_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    type :: distBand_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **spectral peak energy** parameter of the Band spectral model/distribution
    !>  from the corresponding **spectral break energy** \f$\ebreak\f$ and the Band model spectral indices \f$(\alpha, \beta)\f$.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  The break energy of the Band distribution is related to the peak energy via the following relation,
    !>  \f{equation}{
    !>      \epeak = \frac{2 + \alpha}{\alpha - \beta} \ebreak
    !>  \f}
    !>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `ebreak`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `ebreak`, containing the second shape parameter of the distribution.<br>
    !>  \param[in]  ebreak      :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the spectral break energy values.<br>
    !>
    !>  \return
    !>  `epeak`                 :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `ebreak`, containing the peak energy of the distribution.<br>
    !>                              The output value has the same physical units as the input `ebreak`.<br>
    !>
    !>  \interface{getBandEpeak}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: getBandEpeak
    !>
    !>      epeak = getBandEpeak(alpha, beta, ebreak)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [setBandMean](@ref pm_distBand::setBandMean)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>
    !>  \example{getBandEpeak}
    !>  \include{lineno} example/pm_distBand/getBandEpeak/main.F90
    !>  \compilef{getBandEpeak}
    !>  \output{getBandEpeak}
    !>  \include{lineno} example/pm_distBand/getBandEpeak/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{getBandEpeak}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBandEpeak

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBandEpeak_RK5(alpha, beta, ebreak) result(epeak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEpeak_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: epeak
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBandEpeak_RK4(alpha, beta, ebreak) result(epeak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEpeak_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: epeak
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBandEpeak_RK3(alpha, beta, ebreak) result(epeak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEpeak_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: epeak
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBandEpeak_RK2(alpha, beta, ebreak) result(epeak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEpeak_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: epeak
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBandEpeak_RK1(alpha, beta, ebreak) result(epeak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEpeak_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: epeak
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **spectral break energy** parameter of the Band spectral model/distribution
    !>  from the corresponding **spectral peak energy** \f$\epeak\f$ and the Band model spectral indices \f$(\alpha, \beta)\f$.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  The break energy of the Band distribution is related to the peak energy via the following relation,
    !>  \f{equation}{
    !>      \ebreak = \frac{\alpha - \beta}{2 + \alpha} \epeak
    !>  \f}
    !>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `epeak`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `epeak`, containing the second shape parameter of the distribution.<br>
    !>  \param[in]  epeak       :   The input scalar or array of the same shape as other array like arguments,
    !>                              of type `real` of kind \RKALL, containing the spectral peak energy values.<br>
    !>
    !>  \return
    !>  `ebreak`                :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `epeak`, containing the break energy of the distribution.<br>
    !>                              The output value has the same physical units as the input `epeak`.<br>
    !>
    !>  \interface{getBandEbreak}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: getBandEbreak
    !>
    !>      ebreak = getBandEbreak(alpha, beta, epeak)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < epeak` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>
    !>  \example{getBandEbreak}
    !>  \include{lineno} example/pm_distBand/getBandEbreak/main.F90
    !>  \compilef{getBandEbreak}
    !>  \output{getBandEbreak}
    !>  \include{lineno} example/pm_distBand/getBandEbreak/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{getBandEbreak}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBandEbreak

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBandEbreak_RK5(alpha, beta, epeak) result(ebreak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEbreak_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: alpha, beta, epeak
        real(RKG)                   :: ebreak
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBandEbreak_RK4(alpha, beta, epeak) result(ebreak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEbreak_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: alpha, beta, epeak
        real(RKG)                   :: ebreak
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBandEbreak_RK3(alpha, beta, epeak) result(ebreak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEbreak_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: alpha, beta, epeak
        real(RKG)                   :: ebreak
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBandEbreak_RK2(alpha, beta, epeak) result(ebreak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEbreak_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: alpha, beta, epeak
        real(RKG)                   :: ebreak
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBandEbreak_RK1(alpha, beta, epeak) result(ebreak)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandEbreak_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: alpha, beta, epeak
        real(RKG)                   :: ebreak
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **coefficient of continuity** of the Band spectral model/distribution
    !>  from the Band model parameters: the break energy \f$\ebreak\f$ and the Band model spectral indices \f$(\alpha, \beta)\f$.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  The **coefficient of continuity** of the Band model makes the Band distribution **continuously differentiable** and is defined as,
    !>  \f{equation}{
    !>      \large
    !>      \zeta = \left(\frac{\ebreak}{\kev}\right)^{\alpha - \beta} \exp\left(\beta - \alpha\right) ~.
    !>  \f}
    !>  This factor is required for computing the density function of the Band distribution.<br>
    !>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `ebreak`, containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments,
    !>                              of the same type and kind as `ebreak`, containing the second shape parameter of the distribution.<br>
    !>  \param[in]  ebreak      :   The input scalar or array of the same shape as other array like arguments of type `real` of kind \RKALL,
    !>                              containing the normalized (unitless) spectral break energy values: \f$\ebreak = \frac{\ebreak}{\kev}\f$.<br>
    !>
    !>  \return
    !>  `zeta`                  :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `ebreak`, containing the coefficient of continuity of the distribution.<br>
    !>
    !>  \interface{getBandZeta}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: getBandZeta
    !>
    !>      zeta = getBandZeta(alpha, beta, ebreak)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>
    !>  \example{getBandZeta}
    !>  \include{lineno} example/pm_distBand/getBandZeta/main.F90
    !>  \compilef{getBandZeta}
    !>  \output{getBandZeta}
    !>  \include{lineno} example/pm_distBand/getBandZeta/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{getBandZeta}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBandZeta

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBandZeta_RK5(alpha, beta, ebreak) result(zeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandZeta_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: zeta
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBandZeta_RK4(alpha, beta, ebreak) result(zeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandZeta_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: zeta
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBandZeta_RK3(alpha, beta, ebreak) result(zeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandZeta_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: zeta
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBandZeta_RK2(alpha, beta, ebreak) result(zeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandZeta_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: zeta
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBandZeta_RK1(alpha, beta, ebreak) result(zeta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandZeta_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: alpha, beta, ebreak
        real(RKG)                   :: zeta
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **unnormalized** density function (UDF) of the Band spectral model/distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  The **unnormalized unit-less** density function of the Band model can be written as,<br>
    !>  \f{equation}{
    !>      \large
    !>      f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) =
    !>      \begin{cases}
    !>          E^\alpha \exp\left(-\frac{E}{\efold}\right) &,~ \ms{if} & 0 < \ms{lb} \leq E < \ebreak \\
    !>          \zeta E^\beta &,~ \ms{if} & \ebreak \leq E < \ms{ub} < +\infty \\
    !>      \end{cases}
    !>  \f}
    !>  where,
    !>  <ol>
    !>      <li>    \f$E\f$ is the (presumably unitless) energy at which the distribution must be computed,
    !>      <li>    \f$\ebreak\f$ is the (presumably unitless or of the same unit as \f$E\f$) spectral break energy,
    !>      <li>    \f$\efold = \frac{\alpha - \beta}{\ebreak}\f$ is the e-folding energy,
    !>              a measure of the **scale** of the distribution below and above which the distribution
    !>              approaches power-law behavior with exponents \f$\alpha\f$ and \f$\beta\f$ respectively,
    !>      <li>    the factor \f$\eta(\alpha, \beta, \ebreak)\f$ is a normalization constant that properly normalizes the PDF,
    !>      <li>    the factor \f$\zeta = \ebreak^{\alpha - \beta} \exp\left(\beta - \alpha\right)\f$ is
    !>              the **coefficient of continuity** of the distribution that makes the distribution **continuously differentiable**.<br>
    !>      <li>    the constants \f$(\ms{lb}, \ms{ub})\f$ represent the lower and upper bounds of the PDF respectively.<br>
    !>  </ol>
    !>
    !>  \param[in]  energy      :   The input scalar or array of the same shape as other array like arguments of type `real` of kind \RKALL,
    !>                              containing the energy at which the UDF must be computed.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `energy`,
    !>                              containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `energy`,
    !>                              containing the second shape parameter of the distribution.<br>
    !>  \param[in]  ebreak      :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `energy`,
    !>                              containing the spectral break energy values.<br>
    !>  \param[in]  zeta        :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `energy`,
    !>                              containing the containing the [coefficient of continuity](@ref pm_distBand::getBandZeta) of the distribution.<br>
    !>                              (**optional**. default = [getBandZeta(alpha, beta, ebreak)](@ref pm_distBand::getBandZeta). Its presence can expedite the computations.)<br>
    !>                              This means that the output `udf` is computed from the lower tail of the distribution.)
    !>  \param[in]  invEfold    :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `energy`,
    !>                              containing the **inverse** of the e-folding energy of the Band model: \f$\frac{1}{\efold} = \frac{\alpha - \beta}{\ebreak}\f$.<br>
    !>                              (**optional**. default = `(alpha - beta) / ebreak`. Its presence can expedite the computations.)<br>
    !>
    !>  \return
    !>  `udf`                   :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `energy`, containing the distribution UDF.<br>
    !>  \return
    !>
    !>  \interface{getBandUDF}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: getBandUDF
    !>
    !>      udf = getBandUDF(energy, alpha, beta, ebreak, zeta = zeta, invEfold = invEfold)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < energy` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invEfold` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  The condition `ebreak = (alpha - beta) * invEfold` must hold for the corresponding input arguments.<br>
    !>  The condition `zeta = getZeta(alpha, beta, ebreak)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The normalization (and the physical units) of the input `energy` is irrelevant as long as
    !>  the input values `ebreak` and `zeta` are computed in the same physical dimensions and with the same normalizations.<br>
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>
    !>  \example{getBandUDF}
    !>  \include{lineno} example/pm_distBand/getBandUDF/main.F90
    !>  \compilef{getBandUDF}
    !>  \output{getBandUDF}
    !>  \include{lineno} example/pm_distBand/getBandUDF/main.out.F90
    !>  \postproc{getBandUDF}
    !>  \include{lineno} example/pm_distBand/getBandUDF/main.py
    !>  \vis{getBandUDF}
    !>  \image html pm_distBand/getBandUDF/getBandUDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{getBandUDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBandUDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getBandUDF_RK5(energy, alpha, beta, ebreak, zeta, invEfold) result(udf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandUDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: energy, alpha, beta, ebreak
        real(RKG)   , intent(in)    , optional  :: zeta, invEfold
        real(RKG)                               :: udf
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getBandUDF_RK4(energy, alpha, beta, ebreak, zeta, invEfold) result(udf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandUDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: energy, alpha, beta, ebreak
        real(RKG)   , intent(in)    , optional  :: zeta, invEfold
        real(RKG)                               :: udf
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getBandUDF_RK3(energy, alpha, beta, ebreak, zeta, invEfold) result(udf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandUDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: energy, alpha, beta, ebreak
        real(RKG)   , intent(in)    , optional  :: zeta, invEfold
        real(RKG)                               :: udf
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getBandUDF_RK2(energy, alpha, beta, ebreak, zeta, invEfold) result(udf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandUDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: energy, alpha, beta, ebreak
        real(RKG)   , intent(in)    , optional  :: zeta, invEfold
        real(RKG)                               :: udf
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getBandUDF_RK1(energy, alpha, beta, ebreak, zeta, invEfold) result(udf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandUDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: energy, alpha, beta, ebreak
        real(RKG)   , intent(in)    , optional  :: zeta, invEfold
        real(RKG)                               :: udf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **unnormalized** cumulative distribution function (UCDF) of the Band spectral model/distribution.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  The UCDF of the Band model is the integral of its [UDF](@ref pm_distBand::getBandUDF) over a range \f$(\ms{lb}, \ms{ub})\f$ written as,<br>
    !>  \f{equation}{
    !>      \large
    !>      \ms{UCDF}_\ms{BAND} = \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE ~.
    !>  \f}
    !>  where \f$f_{\ms{BAND}}\f$ is the [UDF](@ref pm_distBand::getBandUDF) of the Band distribution.<br>
    !>
    !>  While the integration domain should be ideally \f$[0, +\infty)\f$,
    !>  the arbitrary values of \f$\alpha\f$ and \f$\beta\f$ require finite bounds for the integral to be specified by user to ensure convergence.<br>
    !>
    !>  \param[out] ucdf        :   The output scalar or array of the same shape as any input array-like argument,
    !>                              of type `real` of kind \RKALL as the input argument `ucdf`, containing the distribution UCDF.<br>
    !>  \param[in]  lb          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `ucdf`, representing the lower bound of the Band distribution.<br>
    !>  \param[in]  ub          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                              of the same type and kind as the input argument `ucdf`, representing the upper bound of the Band distribution.<br>
    !>  \param[in]  alpha       :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `ucdf`,
    !>                              containing the first shape parameter of the distribution.<br>
    !>  \param[in]  beta        :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `ucdf`,
    !>                              containing the second shape parameter of the distribution.<br>
    !>  \param[in]  ebreak      :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `ucdf`,
    !>                              containing the **normalized** spectral break energy values: \f$\ebreak = \frac{\ebreak}{100\kev}\f$.<br>
    !>  \param[out] info        :   The output scalar of type `integer` of default kind \IK.<br>
    !>                              On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                              If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm or,
    !>                              to the output error returned by brute force integrator [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                              **An negative output value signifies the lack of convergence and failure to compute the UCDF**.<br>
    !>                              This is likely to happen if the input value for `alpha` or `beta` are too extreme.<br>
    !>
    !>  \interface{setBandUCDF}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: setBandUCDF
    !>
    !>      call setBandUCDF(ucdf, lb, ub, alpha, beta, ebreak, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lb` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invEfold` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  The condition `ebreak = (alpha - beta) * invEfold` must hold for the corresponding input arguments.<br>
    !>  The condition `zeta = getZeta(alpha, beta, ebreak)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The normalization (and the physical units) of the input `energy` is irrelevant as long as
    !>  the input values `ebreak` and `zeta` are computed in the same physical dimensions and with the same normalizations.<br>
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [setBandMean](@ref pm_distBand::setBandMean)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>
    !>  \example{setBandUCDF}
    !>  \include{lineno} example/pm_distBand/setBandUCDF/main.F90
    !>  \compilef{setBandUCDF}
    !>  \output{setBandUCDF}
    !>  \include{lineno} example/pm_distBand/setBandUCDF/main.out.F90
    !>  \postproc{setBandUCDF}
    !>  \include{lineno} example/pm_distBand/setBandUCDF/main.py
    !>  \vis{setBandUCDF}
    !>  \image html pm_distBand/setBandUCDF/setBandUCDF.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{setBandUCDF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBandUCDF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandUCDF_RK5(ucdf, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandUCDF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: ucdf
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandUCDF_RK4(ucdf, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandUCDF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: ucdf
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandUCDF_RK3(ucdf, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandUCDF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: ucdf
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandUCDF_RK2(ucdf, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandUCDF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: ucdf
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandUCDF_RK1(ucdf, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandUCDF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: ucdf
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the mean of the Band distribution for an input set of parameters.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  This generic interface performs one of the following quantities:<br>
    !>  <ol>
    !>      <li>    **The mean of the Band distribution.**<br>
    !>              The mean of the Band distribution with an unknown normalization (amplitude) over an arbitrary range \f$(\ms{lb}, \ms{ub})\f$ is defined by the following integral ratio,<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \mu_E = \frac{ \int_{\ms{lb}}^{\ms{ub}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }{ \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE } ~,
    !>              \f}
    !>              where \f$\mu_E\f$ is the mean of the Band distribution.<br>
    !>              The mean of the Band distribution is particularly important for converting the photon fluence of a Band spectrum to the corresponding energy fluence.<br>
    !>              For example, the energy fluence \f$\sergs\f$ of the Band distribution with the same physical unit as \f$E\f$ can be computed from the corresponding the Band photon fluence \f$\sphot\f$ as,
    !>              \f{equation}{
    !>                  \large
    !>                  \sergs = \sphot \mu_E ~,
    !>              \f}
    !>      <li>    **The generalized mean of the Band distribution.**<br>
    !>              Optionally, this generic interface also computes the above integral with new support \f$(\ms{lbnew}, \ms{ubnew})\f$,
    !>              \f{equation}{
    !>                  \large
    !>                  \mu_E = \frac{ \int_{\ms{lbnew}}^{\ms{ubnew}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }{ \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE } ~,
    !>              \f}
    !>              where \f$\mu_E\f$ is not anymore the common definition of the distribution mean, but a *generalization* of the concept.<br>
    !>              This generalized mean facilitates the computation of the energy or photon fluence over a different range from the range of the original photon or energy fluence.<br>
    !>  </ol>
    !>
    !>  \warning
    !>  The input arguments `lbnew`, `ubnew`, `lb`, `ub`, `ebreak` must all be unit-less
    !>  (without physical dimensions) or all have the same physical units (typically, \f$\kev\f$).<br>
    !>
    !>  \note
    !>  The physical units of the input or output units can be changed via the facilities of module [pm_physUnit](@ref pm_physUnit).<br>
    !>
    !>  \param[out]     mean        :   The output positive scalar or array of the same shape as any input array-like argument, of type `real` of kind \RKALL.<br>
    !>                                  <ol>
    !>                                      <li>    If the input arguments `lbnew` and `ubnew` are missing, then `mean` will contain the mean of the Band distribution on return.<br>
    !>                                      <li>    If the input arguments `lbnew` and `ubnew` are present, then `mean` will contain the ratio of the UCDF of the Band distribution
    !>                                              in the range \f$(\ms{lbnew}, \ms{ubnew})\f$ to the UCDF of the Band distribution in the \f$(\ms{lb}, \ms{ub})\f$ as defined above.<br>
    !>                                  </ol>
    !>  \param[in]      lb          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the output argument `mean`, representing the lower bound of the Band distribution.<br>
    !>  \param[in]      ub          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the output argument `mean`, representing the upper bound of the Band distribution.<br>
    !>  \param[in]      alpha       :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `mean`,
    !>                                  containing the first shape parameter of the distribution.<br>
    !>  \param[in]      beta        :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `mean`,
    !>                                  containing the second shape parameter of the distribution.<br>
    !>  \param[in]      ebreak      :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `mean`,
    !>                                  containing the **normalized** spectral break energy values: \f$\ebreak = \frac{\ebreak}{100\kev}\f$.<br>
    !>  \param[out]     info        :   The output scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm or,
    !>                                  to the output error returned by brute force integrator [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the UCDF**.<br>
    !>                                  This is likely to happen if the input value for `alpha` or `beta` are too extreme.<br>
    !>  \param[in]      lbnew       :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the output argument `mean`, representing the **new** lower bound of the Band distribution.<br>
    !>                                  (**optional**, default = `lb`)
    !>  \param[in]      ubnew       :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the output argument `mean`, representing the **new** upper bound of the Band distribution.<br>
    !>                                  (**optional**, default = `ub`)
    !>
    !>  \interface{setBandMean}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: setBandMean
    !>
    !>      call setBandMean(mean, lb, ub, alpha, beta, ebreak, info)
    !>      call setBandMean(mean, lb, ub, alpha, beta, ebreak, info, lbnew, ubnew)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lb` must hold for the corresponding input arguments.<br>
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lbnew` must hold for the corresponding input arguments.<br>
    !>  The condition `lbnew < ubnew` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < fluence` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lbnew` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ubnew` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < invEfold` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  The condition `ebreak = (alpha - beta) * invEfold` must hold for the corresponding input arguments.<br>
    !>  The condition `zeta = getZeta(alpha, beta, ebreak)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The normalization (and the physical units) of the input `energy` is irrelevant as long as
    !>  the input values `ebreak` and `zeta` are computed in the same physical dimensions and with the same normalizations.<br>
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [setBandMean](@ref pm_distBand::setBandMean)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>
    !>  \example{setBandMean}
    !>  \include{lineno} example/pm_distBand/setBandMean/main.F90
    !>  \compilef{setBandMean}
    !>  \output{setBandMean}
    !>  \include{lineno} example/pm_distBand/setBandMean/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{setBandMean}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBandMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandMeanDef_RK5(mean, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandMeanDef_RK4(mean, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandMeanDef_RK3(mean, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandMeanDef_RK2(mean, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandMeanDef_RK1(mean, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandMeanNew_RK5(mean, lb, ub, alpha, beta, ebreak, lbnew, ubnew, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanNew_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak, lbnew, ubnew
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandMeanNew_RK4(mean, lb, ub, alpha, beta, ebreak, lbnew, ubnew, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanNew_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak, lbnew, ubnew
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandMeanNew_RK3(mean, lb, ub, alpha, beta, ebreak, lbnew, ubnew, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanNew_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak, lbnew, ubnew
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandMeanNew_RK2(mean, lb, ub, alpha, beta, ebreak, lbnew, ubnew, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanNew_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak, lbnew, ubnew
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandMeanNew_RK1(mean, lb, ub, alpha, beta, ebreak, lbnew, ubnew, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandMeanNew_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lb, ub, alpha, beta, ebreak, lbnew, ubnew
        real(RKG)   , intent(out)               :: mean
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the photon integral (the photon fluence in units of photon counts) of the Band model for the given distribution
    !>  parameters from the corresponding energy integral of the distribution (the energy fluence in units of the input break energy).<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  This generic interface computes one of the following quantities:<br>
    !>  <ol>
    !>      <li>    **Conversion of energy fluence to photon fluence.**<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \sphot = \sergs \frac
    !>                  { \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  { \int_{\ms{lb}}^{\ms{ub}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  ~.
    !>              \f}
    !>      <li>    **Conversion of energy fluence to photon fluence in a different energy range.**<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \sphot = \sergs \frac
    !>                  { \int_{\ms{lbnew}}^{\ms{ubnew}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  { \int_{\ms{lb}}^{\ms{ub}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  ~.
    !>              \f}
    !>      <li>    **Conversion of photon fluence to photon fluence in a different energy range.**<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \sphot = \sphot \frac
    !>                  { \int_{\ms{lbnew}}^{\ms{ubnew}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  { \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  ~.
    !>              \f}
    !>  </ol>
    !>
    !>  \warning
    !>  The input arguments `lbnew`, `ubnew`, `energy`, `lb`, `ub`, `ebreak` must all be unit-less
    !>  (without physical dimensions) or all have the same physical units (typically, \f$\kev\f$).<br>
    !>
    !>  \note
    !>  The physical units of the input or output arguments can be changed via the facilities of module [pm_physUnit](@ref pm_physUnit).<br>
    !>
    !>  \param[inout]   photon     :   The input/output positive scalar or array of the same shape as any input array-like argument,
    !>                                  of type `real` of kind \RKALL, containing the photon fluence (UCDF) of the distribution.<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `energy` is present, then `photon` has `intent(out)`.<br>
    !>                                              On output, the value of `photon` represents the photon fluence of the Band model in the original or the newly-specified energy window.<br>
    !>                                      <li>    If the input argument `energy` is missing, then `photon` has `intent(inout)`.<br>
    !>                                              On input, the value of `photon` is used as the original photon fluence to computed the amplitude of the Band model in the old energy window.<br>
    !>                                              On output, the value of `photon` represents the photon fluence of the Band model in the new energy window.<br>
    !>                                  </ol>
    !>  \param[in]      energy      :   The input **positive** scalar or array of the same shape as any input array-like argument, of the same type and kind as the input argument `photon`,
    !>                                  representing the energy fluence (integral) of the Band distribution with the specified input distribution parameters.<br>
    !>                                  (**optional**. If missing, the original fluence is read from the input value of `photon` and is assumed to be the photon fluence.)
    !>  \param[in]      lbnew       :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the **new** lower bound of the Band distribution.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `ubnew` is present and `energy` is missing.)
    !>  \param[in]      ubnew       :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the **new** upper bound of the Band distribution.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `lbnew` is present and `energy` is missing.)
    !>  \param[in]      lb          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the lower bound of the Band distribution.<br>
    !>  \param[in]      ub          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the upper bound of the Band distribution.<br>
    !>  \param[in]      alpha       :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `photon`,
    !>                                  containing the first shape parameter of the distribution.<br>
    !>  \param[in]      beta        :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `photon`,
    !>                                  containing the second shape parameter of the distribution.<br>
    !>  \param[in]      ebreak      :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `photon`,
    !>                                  containing the spectral break energy values.<br>
    !>  \param[out]     info        :   The input scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm or,
    !>                                  to the output error returned by brute force integrator [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the UCDF**.<br>
    !>                                  This is likely to happen if the input value for `alpha` or `beta` are too extreme.<br>
    !>
    !>  \interface{setBandPhoton}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: setBandPhoton
    !>
    !>      call setBandPhoton(photon, energy, lb, ub, alpha, beta, ebreak, info)
    !>      call setBandPhoton(photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
    !>      call setBandPhoton(photon, lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lb` must hold for the corresponding input arguments.<br>
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lbnew` must hold for the corresponding input arguments.<br>
    !>  The condition `lbnew < ubnew` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < photon` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < energy` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lbnew` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ubnew` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The normalization (and the physical units) of the input `energy` is irrelevant as long as
    !>  the other dimensional input values (e.g., `ebreak`) are computed in the same physical dimensions and with the same normalizations.<br>
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [setBandMean](@ref pm_distBand::setBandMean)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>
    !>  \example{setBandPhoton}
    !>  \include{lineno} example/pm_distBand/setBandPhoton/main.F90
    !>  \compilef{setBandPhoton}
    !>  \output{setBandPhoton}
    !>  \include{lineno} example/pm_distBand/setBandPhoton/main.out.F90
    !>  \postproc{setBandPhoton}
    !>  \include{lineno} example/pm_distBand/setBandPhoton/main.py
    !>  \vis{setBandPhoton}
    !>  \image html pm_distBand/setBandPhoton/setBandPhoton.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{setBandPhoton}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBandPhoton

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyOldB_RK5(photon, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyOldB_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyOldB_RK4(photon, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyOldB_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyOldB_RK3(photon, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyOldB_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyOldB_RK2(photon, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyOldB_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyOldB_RK1(photon, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyOldB_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandPhotonFromPhotonNewB_RK5(photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromPhotonNewB_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandPhotonFromPhotonNewB_RK4(photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromPhotonNewB_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandPhotonFromPhotonNewB_RK3(photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromPhotonNewB_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandPhotonFromPhotonNewB_RK2(photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromPhotonNewB_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandPhotonFromPhotonNewB_RK1(photon, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromPhotonNewB_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyNewB_RK5(photon, lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyNewB_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyNewB_RK4(photon, lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyNewB_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyNewB_RK3(photon, lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyNewB_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyNewB_RK2(photon, lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyNewB_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandPhotonFromEnergyNewB_RK1(photon, lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandPhotonFromEnergyNewB_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lbnew, ubnew, energy, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: photon
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the energy integral (the energy fluence in units of the input break energy) of the Band model for the
    !>  given distribution parameters from the corresponding photon integral of the distribution (the photon fluence in units of photon counts).<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distBand](@ref pm_distBand) for more information on the Band distribution.<br>
    !>  This generic interface computes one of the following quantities:<br>
    !>  <ol>
    !>      <li>    **Conversion of photon fluence to energy fluence.**<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \sergs = \sphot \frac
    !>                  { \int_{\ms{lb}}^{\ms{ub}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  { \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  ~.
    !>              \f}
    !>      <li>    **Conversion of photon fluence to energy fluence in a different energy range.**<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \sergs = \sphot \frac
    !>                  { \int_{\ms{lbnew}}^{\ms{ubnew}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  { \int_{\ms{lb}}^{\ms{ub}} f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  ~.
    !>              \f}
    !>      <li>    **Conversion of energy fluence to energy fluence in a different energy range.**<br>
    !>              \f{equation}{
    !>                  \large
    !>                  \sergs = \sergs \frac
    !>                  { \int_{\ms{lbnew}}^{\ms{ubnew}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  { \int_{\ms{lb}}^{\ms{ub}} E ~ f_{\ms{BAND}}(E | \alpha, \beta, \ebreak) dE }
    !>                  ~.
    !>              \f}
    !>  </ol>
    !>
    !>  \warning
    !>  The input arguments `lbnew`, `ubnew`, `energy`, `lb`, `ub`, `ebreak` must all be unit-less
    !>  (without physical dimensions) or all have the same physical units (typically, \f$\kev\f$).<br>
    !>
    !>  \note
    !>  The physical units of the input or output arguments can be changed via the facilities of module [pm_physUnit](@ref pm_physUnit).<br>
    !>
    !>  \param[inout]   energy      :   The input/output positive scalar or array of the same shape as any input array-like argument,
    !>                                  of type `real` of kind \RKALL, containing the energy integral of the distribution (energy fluence).<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `photon` is present, then `energy` has `intent(out)`.<br>
    !>                                              On output, the value of `energy` represents the energy fluence of the Band model in the original or the newly-specified energy window.<br>
    !>                                      <li>    If the input argument `photon` is missing, then `energy` has `intent(inout)`.<br>
    !>                                              On input, the value of `energy` is used to computed the amplitude of the Band model in the old energy window.<br>
    !>                                              On output, the value of `energy` represents the energy fluence of the Band model in the new energy window.<br>
    !>                                  </ol>
    !>  \param[in]      photon      :   The input **positive** scalar or array of the same shape as any input array-like argument, of the same type and kind as the input argument `photon`,
    !>                                  representing the energy integral of the Band distribution (photon counts) with the specified input distribution parameters.<br>
    !>                                  (**optional**. If missing, the original fluence is read from the input value of `energy` and is assumed to have energy units.)
    !>  \param[in]      lbnew       :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the **new** lower bound of the Band distribution.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `ubnew` is present and `energy` is missing.)
    !>  \param[in]      ubnew       :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the **new** upper bound of the Band distribution.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `lbnew` is present and `energy` is missing.)
    !>  \param[in]      lb          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the lower bound of the Band distribution.<br>
    !>  \param[in]      ub          :   The input **positive** scalar or array of the same shape as any input array-like argument,
    !>                                  of the same type and kind as the input argument `photon`, representing the upper bound of the Band distribution.<br>
    !>  \param[in]      alpha       :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `photon`,
    !>                                  containing the first shape parameter of the distribution.<br>
    !>  \param[in]      beta        :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `photon`,
    !>                                  containing the second shape parameter of the distribution.<br>
    !>  \param[in]      ebreak      :   The input scalar or array of the same shape as other array-like arguments of the same type and kind as `photon`,
    !>                                  containing the spectral break energy values.<br>
    !>  \param[out]     info        :   The input scalar of type `integer` of default kind \IK.<br>
    !>                                  On output, it is set to **positive** the number of iterations taken for the series representation of the Gamma function to converge.<br>
    !>                                  If the algorithm fails to converge, then `info` is set to the negative of the number of iterations taken by the algorithm or,
    !>                                  to the negative of the output error returned by brute force integrator [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  **An negative output value signifies the lack of convergence and failure to compute the UCDF**.<br>
    !>                                  This is likely to happen if the input value for `alpha` or `beta` are too extreme.<br>
    !>
    !>  \interface{setBandEnergy}
    !>  \code{.F90}
    !>
    !>      use pm_distBand, only: setBandEnergy
    !>      call setBandEnergy(energy, photon, lb, ub, alpha, beta, ebreak, info)
    !>      call setBandEnergy(energy, lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak, info)
    !>      call setBandEnergy(energy, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lb` must hold for the corresponding input arguments.<br>
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lbnew` must hold for the corresponding input arguments.<br>
    !>  The condition `lbnew < ubnew` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < photon` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < energy` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < lbnew` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ubnew` must hold for the corresponding input arguments.<br>
    !>  The condition `alpha /= -2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < ebreak` must hold for the corresponding input arguments.<br>
    !>  The condition `beta < alpha` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The normalization (and the physical units) of the input `energy` is irrelevant as long as
    !>  the other dimensional input values (e.g., `ebreak`) are computed in the same physical dimensions and with the same normalizations.<br>
    !>
    !>  \see
    !>  [getBandUDF](@ref pm_distBand::getBandUDF)<br>
    !>  [setBandUCDF](@ref pm_distBand::setBandUCDF)<br>
    !>  [setBandMean](@ref pm_distBand::setBandMean)<br>
    !>  [getBandZeta](@ref pm_distBand::getBandZeta)<br>
    !>  [getBandEpeak](@ref pm_distBand::getBandEpeak)<br>
    !>  [getBandEbreak](@ref pm_distBand::getBandEbreak)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>
    !>  \example{setBandEnergy}
    !>  \include{lineno} example/pm_distBand/setBandEnergy/main.F90
    !>  \compilef{setBandEnergy}
    !>  \output{setBandEnergy}
    !>  \include{lineno} example/pm_distBand/setBandEnergy/main.out.F90
    !>  \postproc{setBandEnergy}
    !>  \include{lineno} example/pm_distBand/setBandEnergy/main.py
    !>  \vis{setBandEnergy}
    !>  \image html pm_distBand/setBandEnergy/setBandEnergy.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_distBand](@ref test_pm_distBand)
    !>
    !>  \final{setBandEnergy}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBandEnergy

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonOldB_RK5(energy, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonOldB_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonOldB_RK4(energy, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonOldB_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonOldB_RK3(energy, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonOldB_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonOldB_RK2(energy, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonOldB_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonOldB_RK1(energy, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonOldB_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonNewB_RK5(energy, lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonNewB_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonNewB_RK4(energy, lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonNewB_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonNewB_RK3(energy, lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonNewB_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonNewB_RK2(energy, lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonNewB_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandEnergyFromPhotonNewB_RK1(energy, lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromPhotonNewB_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lbnew, ubnew, photon, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(out)               :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module subroutine setBandEnergyFromEnergyNewB_RK5(energy, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromEnergyNewB_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK4_ENABLED
    impure elemental module subroutine setBandEnergyFromEnergyNewB_RK4(energy, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromEnergyNewB_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK3_ENABLED
    impure elemental module subroutine setBandEnergyFromEnergyNewB_RK3(energy, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromEnergyNewB_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK2_ENABLED
    impure elemental module subroutine setBandEnergyFromEnergyNewB_RK2(energy, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromEnergyNewB_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

#if RK1_ENABLED
    impure elemental module subroutine setBandEnergyFromEnergyNewB_RK1(energy, lbnew, ubnew, lb, ub, alpha, beta, ebreak, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBandEnergyFromEnergyNewB_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: lbnew, ubnew, lb, ub, alpha, beta, ebreak
        real(RKG)   , intent(inout)             :: energy
        integer(IK) , intent(out)               :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Bizarrely and frustratingly, Microsoft Windows Subsystem for Linux Ubuntu with GFortran yields Segmentation faults with internal procedure calls.
! This is so unfortunate. To bypass this issue for now, the following subroutine is implemented as separate submodule
! so that the internal shared parameters can be safely passed as submodule parameters.

!WSL_GFORTRAN_BUG     !>  \brief
!WSL_GFORTRAN_BUG     !> Integrate the Band differential spectrum over the input energy range.
!WSL_GFORTRAN_BUG     !>
!WSL_GFORTRAN_BUG     !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   epk             :   The spectral peak energy in units of [keV].
!WSL_GFORTRAN_BUG     !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   beta            :   The upper spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
!WSL_GFORTRAN_BUG     !> \param[out]  photonFluence   :   The fluence in units of photon counts within the input energy range.
!WSL_GFORTRAN_BUG     !> \param[out]  Err             :   An object of class [err_type](@ref pm_err::err_type) containing error-handling information.
!WSL_GFORTRAN_BUG     subroutine getBandPhoton(lowerLim,upperLim,epk,alpha,beta,tolerance,photonFluence,Err)
!WSL_GFORTRAN_BUG #if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!WSL_GFORTRAN_BUG         !DEC$ ATTRIBUTES DLLEXPORT :: getBandPhoton
!WSL_GFORTRAN_BUG #endif
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG        !use pm_quadRomb, only: getQuadRomb
!WSL_GFORTRAN_BUG !WSL_GFORTRAN_BUG         use QuadPackRK1_pmod, only: qag
!WSL_GFORTRAN_BUG         use pm_err, only: err_type
!WSL_GFORTRAN_BUG         implicit none
!WSL_GFORTRAN_BUG         real(RK)        , intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
!WSL_GFORTRAN_BUG         real(RK)        , intent(out)           :: photonFluence
!WSL_GFORTRAN_BUG         type(err_type)        , intent(out)     :: Err
!WSL_GFORTRAN_BUG         character(*, SK), parameter         :: PROCEDURE_NAME = "@getBandPhoton()"
!WSL_GFORTRAN_BUG         real(RK)                        :: ebrk, alphaPlusTwo
!WSL_GFORTRAN_BUG         real(RK)                        :: thisUpperLim, alphaPlusTwoOverEpk, betaPlusOne
!WSL_GFORTRAN_BUG         real(RK)                        :: alphaMinusBand, coef
!WSL_GFORTRAN_BUG         real(RK)                        :: abserr
!WSL_GFORTRAN_BUG         integer(IK)                     :: neval
!WSL_GFORTRAN_BUG         integer(IK)                     :: ierr
!WSL_GFORTRAN_BUG        !real(RK)                        :: alphaPlusOne, logGammaAlphaPlusOne
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         err%occurred = .false._LK
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         if (lowerLim>=upperLim) then
!WSL_GFORTRAN_BUG             photonFluence = 0._RK
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         ! check if the photon indices are consistent with the mathematical rules
!WSL_GFORTRAN_BUG         if (alpha<beta .or. alpha<-2._RK) then
!WSL_GFORTRAN_BUG             photonFluence = -HUGE_RK
!WSL_GFORTRAN_BUG             err%occurred = .true._LK
!WSL_GFORTRAN_BUG             err%msg = MODULE_NAME//PROCEDURE_NAME//SK_": Error occurred: alpha<beta .or. alpha<-2._RK"
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         ! integrate the spectrum
!WSL_GFORTRAN_BUG         alphaPlusTwo = alpha + 2._RK
!WSL_GFORTRAN_BUG         alphaMinusBand = alpha - beta
!WSL_GFORTRAN_BUG         ebrk = epk*alphaMinusBand/alphaPlusTwo
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         if (lowerLim>ebrk) then
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             ! there is only the high energy component in the photonFluence
!WSL_GFORTRAN_BUG             betaPlusOne = beta + 1._RK
!WSL_GFORTRAN_BUG             coef = ebrk**alphaMinusBand * exp(-alphaMinusBand);
!WSL_GFORTRAN_BUG             photonFluence = coef * ( upperLim**betaPlusOne - lowerLim**betaPlusOne ) / betaPlusOne
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG !#if WSL_ENABLED && CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_START
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         elseif (lowerLim<ebrk) then
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             alphaPlusTwoOverEpk = alphaPlusTwo / epk
!WSL_GFORTRAN_BUG             thisUpperLim = min(upperLim,ebrk)
!WSL_GFORTRAN_BUG             !alphaPlusOne = alpha + 1._RK
!WSL_GFORTRAN_BUG             !if (alpha>-1._RK) then
!WSL_GFORTRAN_BUG             !    logGammaAlphaPlusOne = log_gamma( alphaPlusOne )
!WSL_GFORTRAN_BUG             !    ! use the analytical approach to compute the photonFluence:
!WSL_GFORTRAN_BUG             !    ! https://www.wolframalpha.com/input/?i=integrate+x%5Ea+*+exp(-b*x)
!WSL_GFORTRAN_BUG             !    photonFluence = getGamUpper( exponent = alphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , logGammaExponent = logGammaAlphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , lowerLim = alphaPlusTwoOverEpk * lowerLim &
!WSL_GFORTRAN_BUG             !                            , tolerance = tolerance &
!WSL_GFORTRAN_BUG             !                            ) &
!WSL_GFORTRAN_BUG             !             - getGamUpper( exponent = alphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , logGammaExponent = logGammaAlphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , lowerLim = alphaPlusTwoOverEpk * thisUpperLim &
!WSL_GFORTRAN_BUG             !                            , tolerance = tolerance &
!WSL_GFORTRAN_BUG             !                            )
!WSL_GFORTRAN_BUG             !    photonFluence = photonFluence / alphaPlusTwoOverEpk**alphaPlusOne
!WSL_GFORTRAN_BUG             !else
!WSL_GFORTRAN_BUG                 ! use brute-force integration
!WSL_GFORTRAN_BUG                 call qag( f             = getBandCompLowPhoton  &
!WSL_GFORTRAN_BUG                         , a             = lowerLim              &
!WSL_GFORTRAN_BUG                         , b             = thisUpperLim          &
!WSL_GFORTRAN_BUG                         , epsabs        = 0._RK                 &
!WSL_GFORTRAN_BUG                         , epsrel        = tolerance             &
!WSL_GFORTRAN_BUG                         , key           = 1_IK                  &
!WSL_GFORTRAN_BUG                         , result        = photonFluence         &
!WSL_GFORTRAN_BUG                         , abserr        = abserr                &
!WSL_GFORTRAN_BUG                         , neval         = neval                 &
!WSL_GFORTRAN_BUG                         , ier           = ierr                  &
!WSL_GFORTRAN_BUG                         )
!WSL_GFORTRAN_BUG                 !write(*,*) neval
!WSL_GFORTRAN_BUG                 !call getQuadRomb   ( getFunc       = getBandCompLowPhoton &
!WSL_GFORTRAN_BUG                 !                   , xmin          = lowerLim             &
!WSL_GFORTRAN_BUG                 !                   , xmax          = thisUpperLim         &
!WSL_GFORTRAN_BUG                 !                   , tolerance     = 1.e-7_RK             &
!WSL_GFORTRAN_BUG                 !                   , nRefinement   = 10_IK                &
!WSL_GFORTRAN_BUG                 !                   , photonFluence = photonFluence        &
!WSL_GFORTRAN_BUG                 !                   , ierr          = ierr                 &
!WSL_GFORTRAN_BUG                 !                   )
!WSL_GFORTRAN_BUG                 if (ierr/=0_IK) then
!WSL_GFORTRAN_BUG                     photonFluence = -HUGE_RK
!WSL_GFORTRAN_BUG                     err%occurred = .true._LK
!WSL_GFORTRAN_BUG                     err%stat = ierr
!WSL_GFORTRAN_BUG                     err%msg = MODULE_NAME//PROCEDURE_NAME//SK_": Error occurred at QuadPack routine. Check the error code to identify the root cause."
!WSL_GFORTRAN_BUG                     return
!WSL_GFORTRAN_BUG                 end if
!WSL_GFORTRAN_BUG             !end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             if (upperLim>ebrk) then
!WSL_GFORTRAN_BUG                 ! add the remaining part of the photonFluence from the high-energy component
!WSL_GFORTRAN_BUG                 betaPlusOne = beta + 1._RK
!WSL_GFORTRAN_BUG                 alphaMinusBand = alpha - beta
!WSL_GFORTRAN_BUG                 coef = ebrk**alphaMinusBand * exp(-alphaMinusBand);
!WSL_GFORTRAN_BUG                 photonFluence = photonFluence + coef * ( upperLim**betaPlusOne - ebrk**betaPlusOne ) / betaPlusOne
!WSL_GFORTRAN_BUG                 return
!WSL_GFORTRAN_BUG             end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG     contains
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         pure function getBandCompLowPhoton(energy) result(bandCompLow)
!WSL_GFORTRAN_BUG             implicit none
!WSL_GFORTRAN_BUG             real(RK)        , intent(in)    :: energy
!WSL_GFORTRAN_BUG             real(RK)                :: bandCompLow
!WSL_GFORTRAN_BUG             bandCompLow = energy**alpha * exp(-alphaPlusTwoOverEpk*energy)
!WSL_GFORTRAN_BUG         end function getBandCompLowPhoton
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG !#if WSL_ENABLED && CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_STOP
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG     end subroutine getBandPhoton

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Bizarrely and frustratingly, Microsoft Windows Subsystem for Linux Ubuntu with GFortran yields Segmentation faults with internal procedure calls.
! This is so unfortunate. To bypass this issue for now, the following subroutine is implemented as separate submodule
! so that the internal shared parameters can be safely passed as submodule parameters.

!WSL_GFORTRAN_BUG     !>  \brief
!WSL_GFORTRAN_BUG     !> Integrate the Band differential spectrum over the input energy range in units of the input energy.
!WSL_GFORTRAN_BUG     !>
!WSL_GFORTRAN_BUG     !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   epk             :   The spectral peak energy in units of [keV].
!WSL_GFORTRAN_BUG     !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   beta            :   The upper spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
!WSL_GFORTRAN_BUG     !> \param[out]  energyFluence   :   The fluence in units of the input energy (keV) within the input energy range.
!WSL_GFORTRAN_BUG     !> \param[out]  Err             :   An object of class [err_type](@ref pm_err::err_type) containing error-handling information.
!WSL_GFORTRAN_BUG     subroutine getFluenceEnergy(lowerLim,upperLim,epk,alpha,beta,tolerance,energyFluence,Err)
!WSL_GFORTRAN_BUG #if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!WSL_GFORTRAN_BUG         !DEC$ ATTRIBUTES DLLEXPORT :: getFluenceEnergy
!WSL_GFORTRAN_BUG #endif
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG !WSL_GFORTRAN_BUG         use QuadPackRK1_pmod, only: qag
!WSL_GFORTRAN_BUG         use pm_err, only: err_type
!WSL_GFORTRAN_BUG         implicit none
!WSL_GFORTRAN_BUG         real(RK)        , intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
!WSL_GFORTRAN_BUG         type(err_type)        , intent(out)     :: Err
!WSL_GFORTRAN_BUG         real(RK)        , intent(out)           :: energyFluence
!WSL_GFORTRAN_BUG         character(*, SK), parameter         :: PROCEDURE_NAME = "@getFluenceEnergy()"
!WSL_GFORTRAN_BUG         real(RK)                        :: ebrk, alphaPlusTwo
!WSL_GFORTRAN_BUG         real(RK)                        :: thisUpperLim, alphaPlusTwoOverEpk, betaPlusTwo
!WSL_GFORTRAN_BUG         real(RK)                        :: alphaMinusBand, coef, alphaPlusOne
!WSL_GFORTRAN_BUG         real(RK)                        :: abserr
!WSL_GFORTRAN_BUG         integer(IK)                     :: neval
!WSL_GFORTRAN_BUG         integer(IK)                     :: ierr
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         err%occurred = .false._LK
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         if (lowerLim>=upperLim) then
!WSL_GFORTRAN_BUG             energyFluence = 0._RK
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         ! check if the photon indices are consistent with the mathematical rules
!WSL_GFORTRAN_BUG         if (alpha<beta .or. alpha<-2._RK) then
!WSL_GFORTRAN_BUG             energyFluence = -HUGE_RK
!WSL_GFORTRAN_BUG             err%occurred = .true._LK
!WSL_GFORTRAN_BUG             err%msg = MODULE_NAME//PROCEDURE_NAME//SK_": Error occurred: alpha<beta .or. alpha<-2._RK"
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         ! integrate the spectrum
!WSL_GFORTRAN_BUG         alphaPlusTwo = alpha + 2._RK
!WSL_GFORTRAN_BUG         alphaMinusBand = alpha - beta
!WSL_GFORTRAN_BUG         ebrk = epk*alphaMinusBand/alphaPlusTwo
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         if (lowerLim>ebrk) then
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             ! there is only the high energy component in the energyFluence
!WSL_GFORTRAN_BUG             betaPlusTwo = beta + 2._RK
!WSL_GFORTRAN_BUG             coef = ebrk**alphaMinusBand * exp(-alphaMinusBand);
!WSL_GFORTRAN_BUG             energyFluence = coef * ( upperLim**betaPlusTwo - lowerLim**betaPlusTwo ) / betaPlusTwo
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG !#if WSL_ENABLED && CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_START
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         elseif (lowerLim<ebrk) then
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             alphaPlusTwoOverEpk = alphaPlusTwo / epk
!WSL_GFORTRAN_BUG             thisUpperLim = min(upperLim,ebrk)
!WSL_GFORTRAN_BUG             alphaPlusOne = alpha + 1._RK
!WSL_GFORTRAN_BUG             ! use brute-force integration
!WSL_GFORTRAN_BUG             call qag( f             = getBandCompLowEnergy  &
!WSL_GFORTRAN_BUG                     , a             = lowerLim              &
!WSL_GFORTRAN_BUG                     , b             = thisUpperLim          &
!WSL_GFORTRAN_BUG                     , epsabs        = 0._RK                 &
!WSL_GFORTRAN_BUG                     , epsrel        = tolerance             &
!WSL_GFORTRAN_BUG                     , key           = 1_IK                  &
!WSL_GFORTRAN_BUG                     , result        = energyFluence         &
!WSL_GFORTRAN_BUG                     , abserr        = abserr                &
!WSL_GFORTRAN_BUG                     , neval         = neval                 &
!WSL_GFORTRAN_BUG                     , ier           = ierr                  &
!WSL_GFORTRAN_BUG                     )
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             if (ierr/=0_IK) then
!WSL_GFORTRAN_BUG                 energyFluence = -HUGE_RK
!WSL_GFORTRAN_BUG                 err%occurred = .true._LK
!WSL_GFORTRAN_BUG                 err%stat = ierr
!WSL_GFORTRAN_BUG                 err%msg = MODULE_NAME//PROCEDURE_NAME//SK_": Error occurred at QuadPack routine. Check the error code to identify the root cause."
!WSL_GFORTRAN_BUG                 return
!WSL_GFORTRAN_BUG             end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG             if (upperLim>ebrk) then ! add the remaining part of the energyFluence from the high-energy component
!WSL_GFORTRAN_BUG                 betaPlusTwo = beta + 2._RK
!WSL_GFORTRAN_BUG                 alphaMinusBand = alpha - beta
!WSL_GFORTRAN_BUG                 coef = ebrk**alphaMinusBand * exp(-alphaMinusBand);
!WSL_GFORTRAN_BUG                 energyFluence = energyFluence + coef * ( upperLim**betaPlusTwo - ebrk**betaPlusTwo ) / betaPlusTwo
!WSL_GFORTRAN_BUG                 return
!WSL_GFORTRAN_BUG             end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG     contains
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG         pure function getBandCompLowEnergy(energy) result(bandCompLow)
!WSL_GFORTRAN_BUG             implicit none
!WSL_GFORTRAN_BUG             real(RK)        , intent(in)    :: energy
!WSL_GFORTRAN_BUG             real(RK)                :: bandCompLow
!WSL_GFORTRAN_BUG             bandCompLow = energy**alphaPlusOne * exp(-alphaPlusTwoOverEpk*energy)
!WSL_GFORTRAN_BUG         end function getBandCompLowEnergy
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG !#if WSL_ENABLED && CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_STOP
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG
!WSL_GFORTRAN_BUG     end subroutine getFluenceEnergy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distBand ! LCOV_EXCL_LINE