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
!>  This module contains procedures and generic interfaces for computing the cosmic rates of celestial phenomena.
!>
!>  \details
!>  Specifically, the cosmic rate densities of formations of two celestial objects are currently implemented in this module:
!>  <ol>
!>      <li>    the cosmic Star Formation Rate (SFR) density in the Universe under different modeling scenarios,
!>      <li>    the cosmic Long-duration Gamma-Ray Burst (LGRB) Formation Rate (LGRBFR) density in the Universe as a function of cosmological redshift.
!>  </ol>
!>
!>  Cosmic rate
!>  -----------
!>
!>  The rate of occurrence of a cosmic phenomenon can be readily computed by multiplying its
!>  cosmic rate density with the differential comoving volume element of the Universe at the given redshift.<br>
!>  The following two generic interfaces return the **normalized** differential comoving volume element for a variety of possible input arguments:<br>
!>  <ol>
!>      <li>    [getVolComDiffNormed](@ref pm_cosmology::getVolComDiffNormed)<br>
!>      <li>    [setVolComDiffNormed](@ref pm_cosmology::setVolComDiffNormed)<br>
!>  </ol>
!>  See the generic interfaces above and examples below for computing of the cosmic rates of events from their cosmic rate densities.<br>
!>
!>  \see
!>  [pm_cosmology](@ref pm_cosmology)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_cosmicRate

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_cosmicRate"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Comoving Star Formation Rate (SFR) density for a given redshift \f$z\f$
    !>  based on the estimates of [Hopkins and Beacom et al. (2006)](https://ui.adsabs.harvard.edu/abs/2006ApJ...651..142H/abstract),
    !>  assuming a piecewise power-law fit.<br>
    !>
    !>  \param[in]  logzplus1   :   The input **non-negative** scalar, or array of arbitrary rank, of type `real` of kind \RKALL representing
    !>                              the natural logarithm of the cosmological redshift **plus one**, \f$\log(z+1)\f$,
    !>                              at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output of the same type and kind as the input argument `logzplus1` containing
    !>                              the natural logarithm of the cosmic formation rate density of stars at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityH06
    !>
    !>      logRateDensity = getLogRateDensityH06(logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityL08](@ref pm_cosmicRate::getLogRateDensityL08)<br>
    !>  [getLogRateDensityM14](@ref pm_cosmicRate::getLogRateDensityM14)<br>
    !>  [getLogRateDensityM17](@ref pm_cosmicRate::getLogRateDensityM17)<br>
    !>  [getLogRateDensityF18](@ref pm_cosmicRate::getLogRateDensityF18)<br>
    !>  [getLogRateDensityB10](@ref pm_cosmicRate::getLogRateDensityB10)<br>
    !>  [getLogRateDensityP15](@ref pm_cosmicRate::getLogRateDensityP15)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityH06/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityH06/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityH06/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityH06/getLogRateDensityH06.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityH06/getLogRateDensityH06.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityH06/getLogRateDensityH06.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityH06

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityH06_D0_RK5(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityH06_D0_RK4(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityH06_D0_RK3(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityH06_D0_RK2(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityH06_D0_RK1(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getLogRateDensityH06_D1_RK5(logzplus1) result(logRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)                :: logzplus1
!        real(RKG)                               :: logRateDensity
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getLogRateDensityH06_D1_RK4(logzplus1) result(logRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)                :: logzplus1
!        real(RKG)                               :: logRateDensity
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getLogRateDensityH06_D1_RK3(logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in), contiguous    :: logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getLogRateDensityH06_D1_RK2(logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in), contiguous    :: logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getLogRateDensityH06_D1_RK1(logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in), contiguous    :: logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Comoving Star Formation Rate (SFR) density for a given redshift \f$z\f$
    !>  based on the estimates of [Li (2008)](https://academic-oup-com.ezproxy.uta.edu/mnras/article/388/4/1487/980802),
    !>  assuming a piecewise power-law fit.<br>
    !>
    !>  \param[in]  logzplus1   :   The input **non-negative** scalar, or array of arbitrary rank, of type `real` of kind \RKALL representing
    !>                              the natural logarithm of the cosmological redshift **plus one**, \f$\log(z+1)\f$,
    !>                              at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output of the same type and kind as the input argument `logzplus1` containing
    !>                              the natural logarithm of the cosmic formation rate density of stars at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityL08
    !>
    !>      logRateDensity = getLogRateDensityL08(logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityH06](@ref pm_cosmicRate::getLogRateDensityH06)<br>
    !>  [getLogRateDensityM14](@ref pm_cosmicRate::getLogRateDensityM14)<br>
    !>  [getLogRateDensityM17](@ref pm_cosmicRate::getLogRateDensityM17)<br>
    !>  [getLogRateDensityF18](@ref pm_cosmicRate::getLogRateDensityF18)<br>
    !>  [getLogRateDensityB10](@ref pm_cosmicRate::getLogRateDensityB10)<br>
    !>  [getLogRateDensityP15](@ref pm_cosmicRate::getLogRateDensityP15)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityL08/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityL08/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityL08/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityL08/getLogRateDensityL08.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityL08/getLogRateDensityL08.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityL08/getLogRateDensityL08.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityL08

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityL08_D0_RK5(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityL08_D0_RK4(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityL08_D0_RK3(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityL08_D0_RK2(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityL08_D0_RK1(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getLogRateDensityL08_D1_RK5(logzplus1) result(logRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in)                :: logzplus1
!        real(RKG)                               :: logRateDensity
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getLogRateDensityL08_D1_RK4(logzplus1) result(logRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in)                :: logzplus1
!        real(RKG)                               :: logRateDensity
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getLogRateDensityL08_D1_RK3(logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in), contiguous    :: logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getLogRateDensityL08_D1_RK2(logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in), contiguous    :: logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getLogRateDensityL08_D1_RK1(logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in), contiguous    :: logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Comoving Star Formation Rate (SFR) density for a given redshift \f$z\f$ based on the estimates of
    !>  Eqn. 15 of [Madau 2014: Cosmic Star-Formation History](https://www.annualreviews.org/doi/abs/10.1146/annurev-astro-081811-125615).
    !>
    !>  \details
    !>  \f{equation}{
    !>      \large
    !>      \dot{\rho} = 0.015 \frac{(1+z)^{2.7}} {1 + \big[(1+z)/2.9\big]^{5.6}} ~.
    !>  \f}
    !>
    !>  \param[in]  zplus1      :   The input **non-negative** scalar, or array of the same rank as other array-like arguments, of type `real` of kind \RKALL representing the
    !>                              cosmological redshift **plus one**, \f$(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \param[in]  logzplus1   :   The input **non-negative** scalar, or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the natural logarithm
    !>                              of the redshift **plus one**, \f$\log(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output scalar, or array of the same rank as other array-like arguments, of the same type and kind as the input argument `zplus1`
    !>                              containing the natural logarithm of the cosmic formation rate density of stars at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityM14
    !>
    !>      logRateDensity = getLogRateDensityM14(zplus1, logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityH06](@ref pm_cosmicRate::getLogRateDensityH06)<br>
    !>  [getLogRateDensityL08](@ref pm_cosmicRate::getLogRateDensityL08)<br>
    !>  [getLogRateDensityM17](@ref pm_cosmicRate::getLogRateDensityM17)<br>
    !>  [getLogRateDensityF18](@ref pm_cosmicRate::getLogRateDensityF18)<br>
    !>  [getLogRateDensityB10](@ref pm_cosmicRate::getLogRateDensityB10)<br>
    !>  [getLogRateDensityP15](@ref pm_cosmicRate::getLogRateDensityP15)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityM14/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityM14/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityM14/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityM14/getLogRateDensityM14.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityM14/getLogRateDensityM14.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityM14/getLogRateDensityM14.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityM14

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityM14_D0_RK5(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityM14_D0_RK4(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityM14_D0_RK3(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityM14_D0_RK2(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityM14_D0_RK1(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getLogRateDensityM14_D1_RK5(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getLogRateDensityM14_D1_RK4(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getLogRateDensityM14_D1_RK3(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getLogRateDensityM14_D1_RK2(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getLogRateDensityM14_D1_RK1(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Comoving Star Formation Rate (SFR) density for a given redshift \f$z\f$ based on the estimates of
    !>  Eqn. 15 of [Madau 2014: Cosmic Star-Formation History](https://www.annualreviews.org/doi/abs/10.1146/annurev-astro-081811-125615).
    !>
    !>  \f{equation}{
    !>      \large
    !>      \dot{\rho} = 0.015 \frac{(1+z)^{2.7}} {1 + \big[(1+z)/2.9\big]^{5.6}} ~.
    !>  \f}
    !>
    !>  \param[in]  zplus1      :   The input **non-negative** scalar, or array of the same rank as other array-like arguments, of type `real` of kind \RKALL
    !>                              representing the cosmological redshift **plus one**, \f$(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \param[in]  logzplus1   :   The input **non-negative** scalar, or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the natural logarithm
    !>                              of the redshift **plus one**, \f$\log(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output scalar, or array of the same rank as other array-like arguments, of the same type and kind as the input argument `zplus1`
    !>                              containing the natural logarithm of the cosmic formation rate density of stars at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityM17
    !>
    !>      logRateDensity = getLogRateDensityM17(zplus1, logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityH06](@ref pm_cosmicRate::getLogRateDensityH06)<br>
    !>  [getLogRateDensityL08](@ref pm_cosmicRate::getLogRateDensityL08)<br>
    !>  [getLogRateDensityM14](@ref pm_cosmicRate::getLogRateDensityM14)<br>
    !>  [getLogRateDensityF18](@ref pm_cosmicRate::getLogRateDensityF18)<br>
    !>  [getLogRateDensityB10](@ref pm_cosmicRate::getLogRateDensityB10)<br>
    !>  [getLogRateDensityP15](@ref pm_cosmicRate::getLogRateDensityP15)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityM17/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityM17/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityM17/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityM17/getLogRateDensityM17.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityM17/getLogRateDensityM17.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityM17/getLogRateDensityM17.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityM17

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityM17_D0_RK5(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityM17_D0_RK4(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityM17_D0_RK3(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityM17_D0_RK2(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityM17_D0_RK1(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getLogRateDensityM17_D1_RK5(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getLogRateDensityM17_D1_RK4(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getLogRateDensityM17_D1_RK3(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getLogRateDensityM17_D1_RK2(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getLogRateDensityM17_D1_RK1(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Comoving Star Formation Rate (SFR) density for a given redshift \f$z\f$ based on the estimates of
    !>  [Fermi LAT Collaboration 2018: A gamma-ray determination of the Universe’s star formation history](https://www.science.org/doi/10.1126/science.aat8123).
    !>
    !>  \f{equation}{
    !>      \large
    !>      \dot{\rho} = 0.013 \frac{(1+z)^{2.99}} {1 + \big[(1+z)/2.63\big]^{6.19}} ~.
    !>  \f}
    !>
    !>  \param[in]  zplus1      :   The input **non-negative** scalar, or array of the same rank as other array-like arguments, of type `real` of kind \RKALL
    !>                              representing the cosmological redshift **plus one**, \f$(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \param[in]  logzplus1   :   The input **non-negative** scalar, or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the natural logarithm
    !>                              of the redshift **plus one**, \f$\log(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output scalar, or array of the same rank as other array-like arguments, of the same type and kind as the input argument `zplus1`
    !>                              containing the natural logarithm of the cosmic formation rate density of stars at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityF18
    !>
    !>      logRateDensity = getLogRateDensityF18(zplus1, logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityH06](@ref pm_cosmicRate::getLogRateDensityH06)<br>
    !>  [getLogRateDensityL08](@ref pm_cosmicRate::getLogRateDensityL08)<br>
    !>  [getLogRateDensityM14](@ref pm_cosmicRate::getLogRateDensityM14)<br>
    !>  [getLogRateDensityM17](@ref pm_cosmicRate::getLogRateDensityM17)<br>
    !>  [getLogRateDensityB10](@ref pm_cosmicRate::getLogRateDensityB10)<br>
    !>  [getLogRateDensityP15](@ref pm_cosmicRate::getLogRateDensityP15)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityF18/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityF18/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityF18/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityF18/getLogRateDensityF18.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityF18/getLogRateDensityF18.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityF18/getLogRateDensityF18.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityF18

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityF18_D0_RK5(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityF18_D0_RK4(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityF18_D0_RK3(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityF18_D0_RK2(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityF18_D0_RK1(zplus1, logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1, logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getLogRateDensityF18_D1_RK5(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getLogRateDensityF18_D1_RK4(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getLogRateDensityF18_D1_RK3(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getLogRateDensityF18_D1_RK2(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getLogRateDensityF18_D1_RK1(zplus1, logzplus1) result(LogRateDensity)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)   , intent(in), contiguous    :: zplus1(:), logzplus1(:)
!        real(RKG)                               :: LogRateDensity(size(logzplus1))
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Comoving LGRB Formation Rate (LGFR) density for a given redshift \f$z\f$
    !>  based on the Metallicity-corrected estimates of [Butler et al. (2010): B10](https://iopscience.iop.org/article/10.1088/0004-637X/711/1/495),
    !>  assuming a piecewise power-law LGFR model.
    !>
    !>  \details
    !>  The model of B10 is based on the [Poweto distribution](@ref pm_distPoweto) Star Formation Rate (SFR)
    !>  model proposed by [Hopkins and Beacom (2006)](https://ui.adsabs.harvard.edu/abs/2006ApJ...651..142H/abstract).<br>
    !>
    !>  \f{equation}{
    !>
    !>  \f}
    !>
    !>  \param[in]  logzplus1   :   The input scalar or `contiguous` array of arbitrary rank of type `real` of kind \RKALL containing
    !>                              the natural logarithm of the redshift **plus one**, \f$\log(z+1)\f$,
    !>                              at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output of the same type, kind, and rank as the input argument `logzplus1` containing
    !>                              the natural logarithm of the **unnormalized** cosmic formation rate density of LGRBs at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityB10
    !>
    !>      logRateDensity = getLogRateDensityB10(logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input argument `logzplus1` must be non-negative since negative redshift is cosmologically undefined.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityH06](@ref pm_cosmicRate::getLogRateDensityH06)<br>
    !>  [getLogRateDensityL08](@ref pm_cosmicRate::getLogRateDensityL08)<br>
    !>  [getLogRateDensityM14](@ref pm_cosmicRate::getLogRateDensityM14)<br>
    !>  [getLogRateDensityM17](@ref pm_cosmicRate::getLogRateDensityM17)<br>
    !>  [getLogRateDensityF18](@ref pm_cosmicRate::getLogRateDensityF18)<br>
    !>  [getLogRateDensityP15](@ref pm_cosmicRate::getLogRateDensityP15)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityB10/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityB10/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityB10/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityB10/getLogRateDensityB10.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityB10/getLogRateDensityB10.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityB10/getLogRateDensityB10.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityB10

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityB10_D0_RK5(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityB10_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityB10_D0_RK4(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityB10_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityB10_D0_RK3(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityB10_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityB10_D0_RK2(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityB10_D0_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityB10_D0_RK1(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityB10_D0_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the unnormalized Gamma-Ray Burst Formation Rate (GRBFR) density based on the estimates of Petrosian et al (2015).
    !>
    !>  \param[in]  logzplus1   :   The non-negative input scalar or array of arbitrary rank of type `real` of kind \RKALL containing the natural logarithm of the
    !>                              redshift **plus one**, \f$\log(z+1)\f$, at which the formation rate density must be computed.
    !>
    !>  \return
    !>  `logRateDensity`        :   The output of the same type, kind, and rank as the input argument `logzplus1` containing
    !>                              the natural logarithm of the **unnormalized** formation rate density of LGRBs at the requested redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmicRate, only: getLogRateDensityP15
    !>
    !>      logRateDensity = getLogRateDensityP15(logzplus1)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input argument `logzplus1` must be non-negative since a negative redshift is cosmologically undefined.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogRateDensityB10](@ref pm_cosmicRate::getLogRateDensityB10)<br>
    !>  [getLogRateDensityH06](@ref pm_cosmicRate::getLogRateDensityH06)<br>
    !>  [getLogRateDensityL08](@ref pm_cosmicRate::getLogRateDensityL08)<br>
    !>  [getLogRateDensityM14](@ref pm_cosmicRate::getLogRateDensityM14)<br>
    !>  [getLogRateDensityM17](@ref pm_cosmicRate::getLogRateDensityM17)<br>
    !>  [getLogRateDensityF18](@ref pm_cosmicRate::getLogRateDensityF18)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityP15/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityP15/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmicRate/getLogRateDensityP15/main.py
    !>  \vis
    !>  \image html pm_cosmicRate/getLogRateDensityP15/getLogRateDensityP15.z.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityP15/getLogRateDensityP15.z.sample.png width=700
    !>  \image html pm_cosmicRate/getLogRateDensityP15/getLogRateDensityP15.logzplus1.sample.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmicRate](@ref test_pm_cosmicRate)
    !>
    !>  \todo
    !>  This generic interface can be extended to higher-rank input arrays.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogRateDensityP15

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogRateDensityP15_D0_RK5(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityP15_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogRateDensityP15_D0_RK4(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityP15_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogRateDensityP15_D0_RK3(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityP15_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogRateDensityP15_D0_RK2(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityP15_D0_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogRateDensityP15_D0_RK1(logzplus1) result(logRateDensity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityP15_D0_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)                :: logzplus1
        real(RKG)                               :: logRateDensity
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_cosmicRate