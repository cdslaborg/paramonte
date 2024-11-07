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
!>  This module contains procedures and generic interfaces and constants for cosmological calculations.
!>
!>  \details
!>
!>  Cosmological distances
!>  ----------------------
!>
!>  In addition to common cosmological constants, this module also contains a set of
!>  generic interfaces for computing various measures of cosmological distances.<br>
!>
!>  Specifically, the following measures are currently implemented:<br>
!>  <ul>
!>      <li>    [Universe Age](@ref pm_cosmology::getSizeUnivNormed)<br>
!>      <li>    [Universe Size](@ref pm_cosmology::getSizeUnivNormed)<br>
!>      <li>    [Lookback Time](@ref pm_cosmology::getDisLookbackNormed)<br>
!>      <li>    [Lookback (Light-Travel) Distance](@ref pm_cosmology::getDisLookbackNormed)<br>
!>      <li>    [Angular Diameter Distance](@ref pm_cosmology::getDisAngNormed)<br>
!>      <li>    [Line-of-Sight Comoving Distance](@ref pm_cosmology::getDisComNormed)<br>
!>      <li>    [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed)<br>
!>      <li>    [Luminosity Distance](@ref pm_cosmology::getDisLumNormed)<br>
!>      <li>    [Luminosity / Transverse Comoving Distance Approximation (in LCDM cosmology)](@ref pm_cosmology::getDisComTransNormedWU10)<br>
!>  </ul>
!>
!>  Cosmological volumes
!>  --------------------
!>
!>  The **normalized** comoving volume of the Universe as a function of redshift and cosmological
!>  model parameters can be computed via [getVolComNormed](@ref pm_cosmology::getVolComNormed).<br>
!>
!>  The differential comoving volume of the Universe as a function of redshift and cosmological
!>  model parameters can be computed via [getVolComDiffNormed](@ref pm_cosmology::getVolComDiffNormed).<br>
!>
!>  Hubble parameter
!>  ----------------
!>
!>  The Hubble parameter (the ratio of the rate of change of the scale factor of the universe to the current value of the scale factor)
!>  can be computed using the generic interface [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq) for various cosmological models.<br>
!>
!>  \see
!>  [pm_cosmicRate](@ref pm_cosmicRate)<br>
!>  [Distance measures in cosmology](https://arxiv.org/abs/astro-ph/9905116)<br>
!>
!>  \test
!>  [test_pm_cosmology](@ref test_pm_cosmology)<br>
!>
!>  \todo
!>  \pmed
!>  A visualization comparison of all cosmological distances should be added here.
!>
!>
!>  \test
!>  [test_pm_distExp](@ref test_pm_distExp)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 2:38 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_cosmology

    use pm_kind, only: RKB, IK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_cosmology"

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! The cosmological constants.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    real(RKB) , parameter :: YR2SEC = 31557600._RKB                                             !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** number of seconds in a year, based on the definition of the [light year](@ref LYR2CM). See also [MEAN_SECONDS_PER_YEAR](@ref pm_dateTime::MEAN_SECONDS_PER_YEAR).
    real(RKB) , parameter :: LIGHT_SPEED = 2.99792458e5_RKB                                     !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** speed of light (km/s).
    real(RKB) , parameter :: LYR2CM = 946073047258080000._RKB                                   !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** conversion of one light-year to centimeters.
    real(RKB) , parameter :: PIPC2M = 9.6939420213600000e16_RKB                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** conversion of \f$\pi\f$ Parsecs to meters.
    real(RKB) , parameter :: MPC2CM = PIPC2M * 1.e8_RKB / acos(-1._RKB)                         !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** conversion of one Mega-Parsec (Mpc) to Centimeters.
    real(RKB) , parameter :: MPC2KM = MPC2CM / 1.e5_RKB                                         !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** conversion of one Mega-Parsec (Mpc) to Kilometers.
    real(RKB) , parameter :: MPC2LY = MPC2CM / LYR2CM                                           !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** conversion of one Mega-Parsec (Mpc) to light years `LY`.
    real(RKB) , parameter :: MPC2GLY = MPC2LY / 1.e9_RKB                                        !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** conversion of one Mega-Parsec (Mpc) to Giga light years `GLY`.
    real(RKB) , parameter :: LOG_MPC2CM = log(MPC2CM)                                           !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of [MPC2CM](@ref pm_cosmology::MPC2CM).
    real(RKB) , parameter :: OMEGA_M = 0.31_RKB                                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the current best estimate of the normalized Dark Matter density in the \f$\lambda\f$CDM Cosmology.
    real(RKB) , parameter :: OMEGA_L = 0.69_RKB                                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the current best estimate of the normalized Dark Energy density in the \f$\lambda\f$CDM Cosmology.
    real(RKB) , parameter :: OMEGA_R = 0.0_RKB                                                  !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the current best estimate of the normalized Radiation density in the \f$\lambda\f$CDM Cosmology.
    real(RKB) , parameter :: OMEGA_K = 0.0_RKB                                                  !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the current best estimate of the normalized Curvature density in the \f$\lambda\f$CDM Cosmology.

    real(RKB) , parameter :: HUBBLE_CONST = 67.7_RKB                                            !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the Hubble constant in units of `km/s/Mpc`.
    real(RKB) , parameter :: LOG_LIGHT_SPEED = log(LIGHT_SPEED)                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of the [speed of light](@ref pm_cosmology::LIGHT_SPEED).
    real(RKB) , parameter :: LOG_HUBBLE_CONST = log(HUBBLE_CONST)                               !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of the Hubble constant in units of `km/s/Mpc`.
    real(RKB) , parameter :: INV_HUBBLE_CONST = 1._RKB / HUBBLE_CONST                           !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the inverse of the Hubble constant in units of `Mpc*s/km`.

    real(RKB) , parameter :: HUBBLE_DISTANCE_MPC = LIGHT_SPEED / HUBBLE_CONST                   !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) in units of `Mpc` (defined as the speed of light divided by the Hubble constant).<br>
                                                                                                !!          See also [HUBBLE_DISTANCE_GLY](@ref pm_cosmology::HUBBLE_DISTANCE_GLY).
    real(RKB) , parameter :: HUBBLE_DISTANCE_GLY = HUBBLE_DISTANCE_MPC * MPC2GLY                !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) in units of Giga light years `GLY`.<br>
                                                                                                !!          See also [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    real(RKB) , parameter :: LOG_HUBBLE_DISTANCE_MPC = log(HUBBLE_DISTANCE_MPC)                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) in units of `Mpc`.
    real(RKB) , parameter :: LOG_HUBBLE_DISTANCE_CM = LOG_HUBBLE_DISTANCE_MPC + LOG_MPC2CM      !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) in units of `cm`.

    real(RKB) , parameter :: HUBBLE_VOLUME_MPC3 = HUBBLE_DISTANCE_MPC**3                        !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3) in units of `Mpc^3`.<br>
                                                                                                !!          See also [HUBBLE_VOLUME_GLY3](@ref pm_cosmology::HUBBLE_VOLUME_GLY3).
    real(RKB) , parameter :: HUBBLE_VOLUME_GLY3 = HUBBLE_VOLUME_MPC3 * MPC2GLY**3               !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3) in units of Giga light years cubed `GLY^3`.<br>
                                                                                                !!          See also [HUBBLE_VOLUME_MPC3](@ref pm_cosmology::HUBBLE_VOLUME_MPC3).
    real(RKB) , parameter :: LOG_HUBBLE_VOLUME_MPC3 = log(HUBBLE_VOLUME_MPC3)                   !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3) in units of `Mpc^3`.
    real(RKB) , parameter :: LOG_HUBBLE_VOLUME_CM3 = LOG_HUBBLE_VOLUME_MPC3 + LOG_MPC2CM * 3    !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the natural logarithm of the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3) in units of `cm^3`.

    real(RKB) , parameter :: HUBBLE_TIME_SEC = MPC2KM / HUBBLE_CONST                            !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **exact** Hubble time in units of seconds.<br>
                                                                                                !!          This the time the Universe needed to expand to its present size, assuming that the Hubble constant
                                                                                                !!          has remained unchanged since the Big Bang. It is defined as the reciprocal of the Hubble constant,
                                                                                                !!          \f$\frac{1}{H_0}\f$. See "An Introduction to Modern Cosmology", Liddle 2003, page 57.
    real(RKB) , parameter :: HUBBLE_TIME_GYR = HUBBLE_DISTANCE_GLY                              !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the **approximate** Hubble time in units of Giga-years.<br>
                                                                                                !!          This the time the Universe needed to expand to its present size, assuming that the Hubble constant
                                                                                                !!          has remained unchanged since the Big Bang. It is defined as the reciprocal of the Hubble constant,
                                                                                                !!          \f$\frac{1}{H_0}\f$. See "An Introduction to Modern Cosmology", Liddle 2003, page 57.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>size</b> (or equivalently, the **age**) of the Universe at the desired redshift
    !>  <i>normalized</i> to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) (or equivalently, to [Hubble Time](@ref pm_cosmology::HUBBLE_TIME_GYR)),
    !>  given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the **Universe Size** at a given redshift \f$z\f$
    !>  is simply related to the [Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq) as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_S(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) =
    !>      D_H \int_{z}^{+\infty} \frac{1}{(1+z') E(z'; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz' ~,
    !>  \f}
    !>
    !>  or equivalently, the **Universe Age** to a cosmological object as redshift \f$z\f$
    !>  is simply related to the [Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq) as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      T_A(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) =
    !>      T_H \int_{z}^{+\infty} \frac{1}{(1+z') E(z'; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz' ~,
    !>  \f}
    !>
    !>  where,
    !>      -#  \f$z\f$ is the redshift,
    !>      -#  \f$T_A\f$ is the cosmological Universe Age,
    !>      -#  \f$D_S\f$ is the cosmological Universe Size,
    !>      -#  \f$E(z; \cdots)\f$ is the [dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq),
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$T_H = \frac{1}{H_0}\f$ is the [Hubble Time](@ref pm_cosmology::HUBBLE_TIME_GYR),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>
    !>  The **Universe Size**, also known as the **Light-Travel Size**, is defined as the distance traveled by light since the Big Bang to a given redshift.<br>
    !>  Dividing the **Universe Size** by the [Speed of Light](@ref pm_cosmology::LIGHT_SPEED) yields the **Universe Age**.<br>
    !>  For instance, the radius of the observable universe in this distance measure becomes the age of the universe multiplied
    !>  by the speed of light (1 light year/year), which turns out to be approximately `13.8` billion light years.<br>
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{D_S}{D_H}\f$.<br>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{T_A}{T_H}\f$.<br>
    !>  The default method of computing the Universe Size is numerical integration via [Adaptive Global Gauss-Kronrod 10-21 Quadrature rule](@ref pm_quadPack).<br>
    !>
    !>  \param[in]  zplus1          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                                  at which the **Universe Size** must be computed.<br>
    !>  \param[in]  omegaM          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                                  (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  reltol          :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval           :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err             :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                                  Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                                  See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `sizeUnivNormed`            :   The output scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` containing the cosmological Universe Size (or Universe Age)
    !>                                  at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) (or [Hubble Time](@ref pm_cosmology::HUBBLE_TIME_GYR)).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getSizeUnivNormed
    !>
    !>      sizeUnivNormed = getSizeUnivNormed(zplus1, reltol, neval = neval, err = err)
    !>      sizeUnivNormed = getSizeUnivNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      sizeUnivNormed = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      sizeUnivNormed = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Note that for some cosmological model parameters, the integral in the definition of the Universe Size/Age is divergent.<br>
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLookbackNormed](@ref pm_cosmology::getDisLookbackNormed)<br>
    !>  [getSizeUnivNormed](@ref pm_cosmology::getSizeUnivNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getSizeUnivNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getSizeUnivNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getSizeUnivNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getSizeUnivNormed/getSizeUnivNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getSizeUnivNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getSizeUnivNormedZ_D0_RK5(zplus1, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getSizeUnivNormedZ_D0_RK4(zplus1, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getSizeUnivNormedZ_D0_RK3(zplus1, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getSizeUnivNormedZ_D0_RK2(zplus1, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getSizeUnivNormedZ_D0_RK1(zplus1, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getSizeUnivNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getSizeUnivNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getSizeUnivNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getSizeUnivNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getSizeUnivNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getSizeUnivNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getSizeUnivNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getSizeUnivNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getSizeUnivNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getSizeUnivNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getSizeUnivNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getSizeUnivNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getSizeUnivNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getSizeUnivNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getSizeUnivNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(sizeUnivNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSizeUnivNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: sizeUnivNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>Lookback Distance</b> (or equivalently, the **Lookback Time**) at the desired redshift
    !>  <i>normalized</i> to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) (or equivalently, to [Hubble Time](@ref pm_cosmology::HUBBLE_TIME_GYR)),
    !>  given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the **Lookback Distance** to a
    !>  cosmological object as redshift \f$z\f$ is simply related to the [Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq) as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_T(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) =
    !>      D_H \int_{0}^{z} \frac{1}{(1+z') E(z'; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz' ~,
    !>  \f}
    !>
    !>  or equivalently, the **Lookback Time** to a cosmological object as redshift \f$z\f$
    !>  is simply related to the [Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq) as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      T_L(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) =
    !>      T_H \int_{0}^{z} \frac{1}{(1+z') E(z'; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz' ~,
    !>  \f}
    !>
    !>  where,
    !>      -#  \f$z\f$ is the redshift,
    !>      -#  \f$D_T\f$ is the cosmological Lookback Distance,
    !>      -#  \f$T_L\f$ is the cosmological Lookback Time,
    !>      -#  \f$E(z; \cdots)\f$ is the [dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq),
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$T_H = \frac{1}{H_0}\f$ is the [Hubble Time](@ref pm_cosmology::HUBBLE_TIME_GYR),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>
    !>  The **Lookback Distance**, also known as the **Light-Travel Distance**, is defined as the distance traveled by light from the given redshift to Earth.<br>
    !>  Dividing the **Lookback Distance** by the [Speed of Light](@ref pm_cosmology::LIGHT_SPEED) yields the **Lookback Time**.<br>
    !>  For instance, the radius of the observable universe in this distance measure becomes the age of the universe multiplied
    !>  by the speed of light (1 light year/year), which turns out to be approximately `13.8` billion light years.<br>
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{D_T}{D_H}\f$.<br>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{T_L}{T_H}\f$.<br>
    !>  The default method of computing the Lookback Distance is numerical integration via [Adaptive Global Gauss-Kronrod 10-21 Quadrature rule](@ref pm_quadPack).<br>
    !>
    !>  \param[in]  zplus1          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                                  at which the **Lookback Distance** must be computed.<br>
    !>  \param[in]  omegaM          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                                  (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  reltol          :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval           :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err             :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                                  Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                                  See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `disLookbackNormed`         :   The output scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` containing the cosmological Lookback Distance (or Lookback Time)
    !>                                  at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC) (or [Hubble Time](@ref pm_cosmology::HUBBLE_TIME_GYR)).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getDisLookbackNormed
    !>
    !>      disLookbackNormed = getDisLookbackNormed(zplus1, reltol, neval = neval, err = err)
    !>      disLookbackNormed = getDisLookbackNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      disLookbackNormed = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      disLookbackNormed = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLookbackNormed](@ref pm_cosmology::getDisLookbackNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getDisLookbackNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getDisLookbackNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getDisLookbackNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getDisLookbackNormed/getDisLookbackNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDisLookbackNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLookbackNormedZ_D0_RK5(zplus1, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLookbackNormedZ_D0_RK4(zplus1, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLookbackNormedZ_D0_RK3(zplus1, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLookbackNormedZ_D0_RK2(zplus1, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLookbackNormedZ_D0_RK1(zplus1, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLookbackNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLookbackNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLookbackNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLookbackNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLookbackNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLookbackNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLookbackNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLookbackNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLookbackNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLookbackNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLookbackNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLookbackNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLookbackNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLookbackNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLookbackNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disLookbackNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLookbackNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLookbackNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>line-of-sight Comoving Distance</b> *normalized* to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC), given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the line-of-sight Comoving Distance to a given cosmological
    !>  redshift \f$z\f$ is defined as (e.g., Peebles, 1993),
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_C(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) =
    !>      D_H \int_{0}^{z} \frac{1}{E(z'; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz' ~,
    !>  \f}
    !>
    !>  where
    !>      -#  \f$D_C(z)\f$ is the cosmological line-of-sight Comoving Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$E(z; \cdots)\f$ is the [dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq),
    !>      -#  \f$z\f$ is the redshift,
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{D_C}{D_H}\f$.<br>
    !>  The default method of computing the Comoving Distance is numerical integration via [Adaptive Global Gauss-Kronrod 10-21 Quadrature rule](@ref pm_quadPack).<br>
    !>
    !>  \param[in]  zplus1      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                              at which the **Comoving Distance** must be computed.<br>
    !>  \param[in]  omegaM      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                              (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  reltol      :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                              A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval       :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                              (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err         :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                              Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                              See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `disComNormed`          :   The output scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` containing the cosmological line-of-sight Comoving Distance
    !>                              at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface{getDisComNormed}
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getDisComNormed
    !>
    !>      disComNormed = getDisComNormed(zplus1, reltol, neval = neval, err = err)
    !>      disComNormed = getDisComNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      disComNormed = getDisComNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      disComNormed = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example{getDisComNormed}
    !>  \include{lineno} example/pm_cosmology/getDisComNormed/main.F90
    !>  \compilef{getDisComNormed}
    !>  \output{getDisComNormed}
    !>  \include{lineno} example/pm_cosmology/getDisComNormed/main.out.F90
    !>  \postproc{getDisComNormed}
    !>  \include{lineno} example/pm_cosmology/getDisComNormed/main.py
    !>  \vis{getDisComNormed}
    !>  \image html pm_cosmology/getDisComNormed/getDisComNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final{getDisComNormed}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDisComNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComNormedZ_D0_RK5(zplus1, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComNormedZ_D0_RK4(zplus1, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComNormedZ_D0_RK3(zplus1, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComNormedZ_D0_RK2(zplus1, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComNormedZ_D0_RK1(zplus1, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>Transverse Comoving Distance</b> <i>normalized</i> to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC), given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the Comoving Distance between two
    !>  events at the same redshift or distance but separated on the sky by some angle \f$\delta\theta\f$ is \f$D_M\delta\theta\f$ where \f$D_M\f$ is the
    !>  Transverse Comoving Distance and is simply related to the [line-of-sight Comoving Distance](@ref pm_cosmology::getDisComNormed) as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_M(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) =
    !>      \begin{cases}
    !>          \frac{D_H}{\sqrt{\Omega_K}} \sinh \bigg( \frac{\sqrt{\Omega_K}}{D_H} D_C \bigg) &, \Omega_K > 0 ~, \\
    !>          D_C &, \Omega_K = 0 ~, \\
    !>          \frac{D_H}{\sqrt{-\Omega_K}} \sin \bigg( \frac{\sqrt{-\Omega_K}}{D_H} D_C \bigg) &, \Omega_K < 0 ~,
    !>      \end{cases}
    !>  \f}
    !>
    !>  where
    !>      -#  \f$D_C(z)\f$ is the cosmological line-of-sight Comoving Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$z\f$ is the redshift,
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{D_M}{D_H}\f$.<br>
    !>  The default method of computing the Transverse Comoving Distance is the same as that of the [line-of-sight Comoving Distance](@ref pm_cosmology::getDisComNormed).<br>
    !>
    !>  \param[in]  zplus1          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                                  at which the **Comoving Distance** must be computed.<br>
    !>  \param[in]  omegaM          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                                  (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  sqrtAbsOmegaK   :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the square root of the normalized curvature density of the universe `sqrt(omegaK)`.<br>
    !>                                  This argument is requested as an input to avoid the redundant costly `sqrt(OmegaK)` operation within the algorithm
    !>                                  when the procedure is to be called repeatedly with the same `omegaK`.<br>
    !>                                  (**optional**, it must be present <b>if and only if</b> `omegaK` is also present.)
    !>  \param[in]  reltol          :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval           :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err             :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                                  Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                                  See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `disComTransNormed`         :   The output scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` containing the cosmological Transverse Comoving Distance
    !>                                  at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getDisComTransNormed
    !>
    !>      disComTransNormed = getDisComTransNormed(zplus1, reltol, neval = neval, err = err)
    !>      disComTransNormed = getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      disComTransNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      disComTransNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getDisComNormed](@ref pm_cosmology::getDisComNormed) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getDisComTransNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getDisComTransNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getDisComTransNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getDisComTransNormed/getDisComTransNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDisComTransNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComTransNormedZ_D0_RK5(zplus1, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComTransNormedZ_D0_RK4(zplus1, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComTransNormedZ_D0_RK3(zplus1, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComTransNormedZ_D0_RK2(zplus1, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComTransNormedZ_D0_RK1(zplus1, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComTransNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComTransNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComTransNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComTransNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComTransNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComTransNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComTransNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComTransNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComTransNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComTransNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisComTransNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisComTransNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisComTransNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisComTransNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisComTransNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disComTransNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disComTransNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>Angular Diameter Distance</b> <i>normalized</i> to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC), given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the **Angular Diameter Distance** to a
    !>  cosmological object as redshift \f$z\f$ is simply related to the [Transverse Comoving Distance)](@ref pm_cosmology::getDisComTransNormed) as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_A(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) = \frac{D_M(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)}{z + 1}
    !>  \f}
    !>
    !>  where
    !>      -#  \f$D_A(z)\f$ is the cosmological Angular Diameter Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$D_M(z)\f$ is the cosmological Transverse Comoving Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$z\f$ is the redshift,
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{D_A}{D_H}\f$.<br>
    !>  The default method of computing the Angular Diameter Distance is the same as that of the [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed).<br>
    !>
    !>  \param[in]  zplus1          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                                  at which the **Angular Diameter Distance** must be computed.<br>
    !>  \param[in]  omegaM          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                                  (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  sqrtAbsOmegaK   :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the square root of the normalized curvature density of the universe `sqrt(omegaK)`.<br>
    !>                                  This argument is requested as an input to avoid the redundant costly `sqrt(OmegaK)` operation within the algorithm
    !>                                  when the procedure is to be called repeatedly with the same `omegaK`.<br>
    !>                                  (**optional**, it must be present <b>if and only if</b> `omegaK` is also present.)
    !>  \param[in]  reltol          :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval           :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err             :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                                  Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                                  See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `disAngNormed`              :   The output scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` containing the cosmological Angular Diameter Distance
    !>                                  at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getDisAngNormed
    !>
    !>      disAngNormed = getDisAngNormed(zplus1, reltol, neval = neval, err = err)
    !>      disAngNormed = getDisAngNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      disAngNormed = getDisAngNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      disAngNormed = getDisAngNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getDisAngNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getDisAngNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getDisAngNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getDisAngNormed/getDisAngNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDisAngNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisAngNormedZ_D0_RK5(zplus1, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisAngNormedZ_D0_RK4(zplus1, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisAngNormedZ_D0_RK3(zplus1, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisAngNormedZ_D0_RK2(zplus1, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisAngNormedZ_D0_RK1(zplus1, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisAngNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisAngNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisAngNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisAngNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisAngNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisAngNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisAngNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisAngNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisAngNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisAngNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisAngNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisAngNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisAngNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisAngNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisAngNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disAngNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisAngNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disAngNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>Luminosity Distance</b> <i>normalized</i> to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC), given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the **Luminosity Distance** to a
    !>  cosmological object as redshift \f$z\f$ is simply related to the [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed) as,
    !>
    !>  \f{eqnarray}{
    !>      \large
    !>      D_L(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)
    !>      &=& (z + 1) ~ D_M(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) ~, \\
    !>      &=& (z + 1)^2 ~ D_A(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) ~,
    !>  \f}
    !>
    !>  where,
    !>      -#  \f$D_L(z; \cdots)\f$ is the cosmological Luminosity Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$D_M(z; \cdots)\f$ is the cosmological Transverse Comoving Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$D_A(z; \cdots)\f$ is the cosmological Angular Diameter Distance in units of `Mpc` as a function of redshift,
    !>      -#  \f$z\f$ is the redshift,
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{D_L}{D_H}\f$.<br>
    !>  The default method of computing the Luminosity Distance is the same as that of the [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed).<br>
    !>
    !>  \param[in]  zplus1          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                                  at which the **Luminosity Distance** must be computed.<br>
    !>  \param[in]  omegaM          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                                  (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK          :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                                  (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  sqrtAbsOmegaK   :   The input scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` representing the square root of the normalized curvature density of the universe `sqrt(omegaK)`.<br>
    !>                                  This argument is requested as an input to avoid the redundant costly `sqrt(OmegaK)` operation within the algorithm
    !>                                  when the procedure is to be called repeatedly with the same `omegaK`.<br>
    !>                                  (**optional**, it must be present <b>if and only if</b> `omegaK` is also present.)
    !>  \param[in]  reltol          :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval           :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                                  (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err             :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                                  Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                                  See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `disLumNormed`              :   The output scalar or array of the same rank as other array-like arguments,
    !>                                  of the same type and kind as the input argument `zplus1` containing the cosmological Luminosity Distance
    !>                                  at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getDisLumNormed
    !>
    !>      disLumNormed = getDisLumNormed(zplus1, reltol, neval = neval, err = err)
    !>      disLumNormed = getDisLumNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      disLumNormed = getDisLumNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      disLumNormed = getDisLumNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getDisLumNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getDisLumNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getDisLumNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getDisLumNormed/getDisLumNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDisLumNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLumNormedZ_D0_RK5(zplus1, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLumNormedZ_D0_RK4(zplus1, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLumNormedZ_D0_RK3(zplus1, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLumNormedZ_D0_RK2(zplus1, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLumNormedZ_D0_RK1(zplus1, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLumNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLumNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLumNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLumNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLumNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLumNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLumNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLumNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLumNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLumNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getDisLumNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getDisLumNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getDisLumNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getDisLumNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getDisLumNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(disLumNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisLumNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: disLumNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **approximate** to the cosmological <b>Transverse Comoving Distance</b> <i>normalized</i> to [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>  given the user-specified cosmological parameters with negligible Curvature (\f$\Omega_K\f$) and Radiation (\f$\Omega_R\f$) density parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R = 0, \Omega_K = 0)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the Transverse Comoving Distance to a given cosmological
    !>  redshift \f$z\f$ is defined as (e.g., Peebles, 1993),
    !>
    !>  \f{equation}{
    !>      \large
    !>      D_M(z; \Omega_M, \Omega_\Lambda) = \int_{0}^{z} \frac{D_H}{\sqrt{(1+z')^3\Omega_M + \Omega_\Lambda}} ~ dz' ~.
    !>  \f}
    !>
    !>  where
    !>      -#  \f$D_M(z; \cdots)\f$ is the cosmological Transverse Comoving Distance as a function of redshift,
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>      -#  \f$z\f$ is the redshift,
    !>
    !>  The approximation used in this implementation is based on the work of
    !>  [Wickramasinghe and Ukwatta (2010), WU10](https://academic.oup.com/mnras/article/406/1/548/1076223?login=true).<br>
    !>  The method of WU10 claims a relative numerical uncertainty of \f$<0.0022\f$ in the computed Luminosity Distance (hence for the Transverse Comoving Distance)
    !>  for \f$(\Omega_M, \Omega_\Lambda) = (0.7, 0.3)\f$.<br>
    !>
    !>  \warning
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{\tilde{D}_M}{D_H}\f$, where \f$\tilde{D}_M\f$ is the computed approximation to the true \f$D_M\f$.<br>
    !>
    !>  \note
    !>  To compute the approximate normalized luminosity distance, simply multiply the output of this routine with the factor `(z + 1)` where `z` stands for redshift.<br>
    !>
    !>  Approximation error
    !>  -------------------
    !>
    !>  The following plots depict the achieved **relative percentage errors** in the computation of the Transverse Comoving Distance for different values of \f$\Omega_\Lambda\f$ and redshift \f$z\f$.<br>
    !>
    !>  \image html pm_cosmology@getDisComTransNormedWU10.png width=900
    !>
    !>  \warning
    !>
    !>  \param[in]  zplus1              :   The input scalar or array of arbitrary rank of type `real` of kind \RKALL containing
    !>                                      the redshift **plus one**, \f$(z+1)\f$, at which the Transverse Comoving Distance must be computed.
    !>  \param[in]  omegaM              :   The input scalar or array of the same rank as other array-like arguments,
    !>                                      of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                                      (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL              :   The input scalar or array of the same rank as other array-like arguments,
    !>                                      of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                                      (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaM` is also present.)
    !>
    !>  \return
    !>  `disComTransNormedWU10`         :   The output of the same type, kind, and rank as the input argument `zplus1` containing the approximate to the Transverse Comoving Distance
    !>                                      at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getDisComTransNormedWU10
    !>
    !>      disComTransNormedWU10 = getDisComTransNormedWU10(zplus1)
    !>      disComTransNormedWU10 = getDisComTransNormedWU10(zplus1, omegaM, omegaL)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Beware that the actual precision of the output of the procedures under this generic interface
    !>  is upper bounded by the approximation algorithm at around relative accuracy of `0.002`.<br>
    !>  This is the case regardless of the `real` kind and precision requested.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getDisComTransNormedWU10/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getDisComTransNormedWU10/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getDisComTransNormedWU10/main.py
    !>  \vis
    !>  \image html pm_cosmology/getDisComTransNormedWU10/getDisComTransNormedWU10.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \todo
    !>  \phigh
    !>  The origins of the glitch in the output of this algorithm at high redshifts (visible in the plots) must be investigated and corrected.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDisComTransNormedWU10

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getDisComTransNormedWU10Z_D0_RK5(zplus1) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10Z_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getDisComTransNormedWU10Z_D0_RK4(zplus1) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10Z_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getDisComTransNormedWU10Z_D0_RK3(zplus1) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10Z_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getDisComTransNormedWU10Z_D0_RK2(zplus1) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10Z_D0_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getDisComTransNormedWU10Z_D0_RK1(zplus1) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10Z_D0_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getDisComTransNormedWU10ZML_D0_RK5(zplus1, omegaM, omegaL) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10ZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getDisComTransNormedWU10ZML_D0_RK4(zplus1, omegaM, omegaL) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10ZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getDisComTransNormedWU10ZML_D0_RK3(zplus1, omegaM, omegaL) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10ZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getDisComTransNormedWU10ZML_D0_RK2(zplus1, omegaM, omegaL) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10ZML_D0_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getDisComTransNormedWU10ZML_D0_RK1(zplus1, omegaM, omegaL) result(disComTransNormedWU10)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisComTransNormedWU10ZML_D0_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)                                   :: disComTransNormedWU10
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **square** of the **dimensionless** Hubble Parameter \f$E(z)^2 = \big(\frac{H(z)}{H_0}\big)^2\f$
    !>  for the default or the specified cosmological parameters.<br>
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the **dimensionless** Hubble Parameter at a given cosmological
    !>  redshift \f$z\f$ is defined as (Peebles, 1993, Principles of Physical Cosmology, pp 310-321),
    !>  \f{equation}{
    !>      \large
    !>      E(z) = \frac{H(z)}{H_0} = \sqrt{\Omega_R(1+z)^4 + \Omega_M(1+z)^3 + \Omega_K(1+z)^2 + \Omega_\Lambda} ~,
    !>  \f}
    !>  where \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).<br>
    !>  Note that \f$E(Z)\f$ is the time derivative of the logarithm of the scale factor \f$a(t)\f$ of the Universe.<br>
    !>
    !>  \param[in]  zplus1      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                              at which the square of the **dimensionless** Hubble Parameter must be computed.<br>
    !>  \param[in]  omegaM      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                              (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>
    !>  \return
    !>  `hubbleParamNormedSq`   :   The output scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` containing the
    !>                              square of the dimensionless Hubble Parameter at the desired redshift.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getHubbleParamNormedSq
    !>
    !>      hubbleParamNormedSq = getHubbleParamNormedSq(zplus1)
    !>      hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL)
    !>      hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)
    !>      hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `omegaM + omegaL + omegaR + omegaK = 1` must hold in all circumstances.<br>
    !>  The condition `1 <= zplus1` must hold in all circumstances.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  \elemental
    !>
    !>  \remark
    !>  The complexity of the output value (instead of simply outputting the Hubble Parameter) is
    !>  to ensure the highest computational efficiency of this and other dependent algorithms by minimizing the potentially redundant calculations.<br>
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Note that `omegaK` (i.e., the **spatial curvature density** of the Universe) is defined as `omegaK = 1 - omegaM - omegaL - omegaR`.<br>
    !>  However, it is requested explicitly in this
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getHubbleParamNormedSq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getHubbleParamNormedSq/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getHubbleParamNormedSq/main.py
    !>  \vis
    !>  \image html pm_cosmology/getHubbleParamNormedSq/getHubbleParamNormedSq.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday 5:43 PM, December 25, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface getHubbleParamNormedSq

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getHubbleParamNormedSqZ_D0_RK5(zplus1) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getHubbleParamNormedSqZ_D0_RK4(zplus1) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getHubbleParamNormedSqZ_D0_RK3(zplus1) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getHubbleParamNormedSqZ_D0_RK2(zplus1) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getHubbleParamNormedSqZ_D0_RK1(zplus1) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getHubbleParamNormedSqZML_D0_RK5(zplus1, omegaM, omegaL) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getHubbleParamNormedSqZML_D0_RK4(zplus1, omegaM, omegaL) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getHubbleParamNormedSqZML_D0_RK3(zplus1, omegaM, omegaL) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getHubbleParamNormedSqZML_D0_RK2(zplus1, omegaM, omegaL) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getHubbleParamNormedSqZML_D0_RK1(zplus1, omegaM, omegaL) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getHubbleParamNormedSqZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHubbleParamNormedSqZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKG)                               :: hubbleParamNormedSq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>Comoving Volume Element</b> per **unit solid angle** of the sky (i.e., `1` [Steradian](https://en.wikipedia.org/wiki/Steradian)
    !>  *normalized* to [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3), given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the Comoving Volume Element \f$dV_C\f$
    !>  at a given cosmological redshift \f$z\f$ is defined as (e.g., Peebles, 1993),
    !>
    !>  \f{eqnarray}{
    !>      \large
    !>      dV_C(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)
    !>      &=& D_H \frac{\big[D_M(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\big]^2}{E(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz ~ d\Omega ~, \\
    !>      &=& D_H \frac{\big[D_A(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\big]^2}{E(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} (1+z)^2 ~ dz ~ d\Omega ~, \\
    !>      &=& D_H \frac{\big[D_L(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\big]^2}{E(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} \frac{1}{(1+z)^2} ~ dz ~ d\Omega ~,
    !>  \f}
    !>
    !>  where
    !>      -#  \f$dV_C(z; \cdots)\f$ is the cosmological Comoving Volume Element at redshift \f$z\f$,,
    !>      -#  \f$D_M(z; \cdots)\f$ is the cosmological [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed) as a function of redshift,
    !>      -#  \f$D_A(z; \cdots)\f$ is the cosmological [Angular Diameter Distance](@ref pm_cosmology::getDisAngNormed) as a function of redshift,
    !>      -#  \f$D_L(z; \cdots)\f$ is the cosmological [Luminosity Distance](@ref pm_cosmology::getDisLumNormed) as a function of redshift,
    !>      -#  \f$E(z; \cdots)\f$ is the [dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq),
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>      -#  \f$z\f$ is the redshift,
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{dV_C}{V_H}\f$, that is, **normalized** to the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3).<br>
    !>  To obtain the **full-sky (normalized) Comoving Volume Element**, multiply the output of the procedures under this generic interface by \f$4\pi\f$, that is, `4 * acos(-1.)`.<br>
    !>  The default method of computing the Comoving Volume Element is the same as that of [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed).<br>
    !>
    !>  \param[in]  zplus1      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of type `real` of kind \RKALL representing the redshift **plus one**, \f$\log(z+1)\f$,
    !>                              at which the **Comoving Volume Element** must be computed.<br>
    !>  \param[in]  omegaM      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized matter density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_M](@ref pm_cosmology::OMEGA_M). It must be present if `omegaL` (\f$\Omega_\Lambda\f$) is also present.)
    !>  \param[in]  omegaL      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized Dark Energy density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_L](@ref pm_cosmology::OMEGA_L). It must be present if `omegaR` is also present.)
    !>  \param[in]  omegaR      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized radiation density in the universe.<br>
    !>                              (**optional**, default = [OMEGA_R](@ref pm_cosmology::OMEGA_R). It must be present if `omegaK` is also present.)
    !>  \param[in]  omegaK      :   The input scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` representing the normalized curvature density of the universe.<br>
    !>                              (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K))
    !>  \param[in]  reltol      :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                              A reasonable recommended value is `reltol = sqrt(epsilon(real(0, kind(zplus1))))`.<br>
    !>  \param[out] neval       :   See the description of the corresponding argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>                              (**optional**. If missing, the number of function evaluations will not be returned.)<br>
    !>  \param[out] err         :   The output scalar of type `integer` of default kind \IK, that is set to zero if the integration converges without any errors.<br>
    !>                              Otherwise, a non-zero value of `err` indicates the occurrence of an error of varying severities.<br>
    !>                              See the description of the corresponding output argument in the documentation of [getQuadErr](@ref pm_quadPack::getQuadErr).<br>
    !>
    !>  \return
    !>  `volComDiffNormed`     :   The output scalar or array of the same rank as other array-like arguments,
    !>                              of the same type and kind as the input argument `zplus1` containing the cosmological Comoving Volume Element
    !>                              at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getVolComDiffNormed
    !>
    !>      volComDiffNormed = getVolComDiffNormed(zplus1, reltol, neval = neval, err = err)
    !>      volComDiffNormed = getVolComDiffNormed(zplus1, omegaM, omegaL, reltol, neval = neval, err = err)
    !>      volComDiffNormed = getVolComDiffNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval = neval, err = err)
    !>      volComDiffNormed = getVolComDiffNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval = neval, err = err)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions listed in the documentation of [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq) must hold for this generic interface.<br>
    !>  \vericons
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  This generic interface performs a numerical integration to compute the output.<br>
    !>  This can be costly for repeated calculations, particularly when the Transverse Comoving distance at the given redshift is known a priori.<br>
    !>  See [setVolComDiffNormed](@ref pm_cosmology::setVolComDiffNormed) for an equivalent, less flexible, but more performant interface that takes advantage of a priori knowledge.<br>
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to the \f$\Lambda\f$CDM Universe with the latest experimental parameter inferences.<br>
    !>
    !>  \note
    !>  Setting `omegaM = 1` and `omegaL = 0` corresponds to the
    !>  [Einstein–de Sitter model of the universe](https://en.wikipedia.org/wiki/Einstein%E2%80%93de_Sitter_universe)
    !>  proposed by Albert Einstein and Willem de Sitter in 1932.<br>
    !>
    !>  \see
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getVolComDiffNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getVolComDiffNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getVolComDiffNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getVolComDiffNormed/getVolComDiffNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getVolComDiffNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getVolComDiffNormedZ_D0_RK5(zplus1, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZ_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getVolComDiffNormedZ_D0_RK4(zplus1, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZ_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getVolComDiffNormedZ_D0_RK3(zplus1, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZ_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getVolComDiffNormedZ_D0_RK2(zplus1, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZ_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getVolComDiffNormedZ_D0_RK1(zplus1, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZ_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getVolComDiffNormedZML_D0_RK5(zplus1, omegaM, omegaL, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZML_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getVolComDiffNormedZML_D0_RK4(zplus1, omegaM, omegaL, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZML_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getVolComDiffNormedZML_D0_RK3(zplus1, omegaM, omegaL, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZML_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getVolComDiffNormedZML_D0_RK2(zplus1, omegaM, omegaL, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZML_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getVolComDiffNormedZML_D0_RK1(zplus1, omegaM, omegaL, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZML_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getVolComDiffNormedZMLR_D0_RK5(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLR_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getVolComDiffNormedZMLR_D0_RK4(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLR_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getVolComDiffNormedZMLR_D0_RK3(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLR_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getVolComDiffNormedZMLR_D0_RK2(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLR_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getVolComDiffNormedZMLR_D0_RK1(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLR_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getVolComDiffNormedZMLRK_D0_RK5(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLRK_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getVolComDiffNormedZMLRK_D0_RK4(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLRK_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getVolComDiffNormedZMLRK_D0_RK3(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLRK_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getVolComDiffNormedZMLRK_D0_RK2(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLRK_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getVolComDiffNormedZMLRK_D0_RK1(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) result(volComDiffNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComDiffNormedZMLRK_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK
        real(RKG)   , intent(in)                    :: reltol
        integer(IK) , intent(out)   , optional      :: neval
        integer(IK) , intent(out)   , optional      :: err
        real(RKG)                                   :: volComDiffNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the cosmological <b>Comoving Volume Differential (Element)</b> per **unit solid angle** of the sky (i.e., `1` [Steradian](https://en.wikipedia.org/wiki/Steradian)
    !>  *normalized* to [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3), given the user-specified cosmological parameters.
    !>
    !>  \details
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the Comoving Volume Element \f$dV_C\f$
    !>  at a given cosmological redshift \f$z\f$ is defined as (e.g., Peebles, 1993),
    !>
    !>  \f{eqnarray}{
    !>      \large
    !>      dV_C(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)
    !>      &=& D_H \frac{\big[D_M(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\big]^2}{E(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} ~ dz ~ d\Omega ~, \\
    !>      &=& D_H \frac{\big[D_A(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\big]^2}{E(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} (1+z)^2 ~ dz ~ d\Omega ~, \\
    !>      &=& D_H \frac{\big[D_L(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\big]^2}{E(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)} \frac{1}{(1+z)^2} ~ dz ~ d\Omega ~,
    !>  \f}
    !>
    !>  where
    !>      -#  \f$dV_C(z; \cdots)\f$ is the cosmological Comoving Volume Element at redshift \f$z\f$,,
    !>      -#  \f$D_M(z; \cdots)\f$ is the cosmological [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed) as a subroutine of redshift,
    !>      -#  \f$D_A(z; \cdots)\f$ is the cosmological [Angular Diameter Distance](@ref pm_cosmology::getDisAngNormed) as a subroutine of redshift,
    !>      -#  \f$D_L(z; \cdots)\f$ is the cosmological [Luminosity Distance](@ref pm_cosmology::getDisLumNormed) as a subroutine of redshift,
    !>      -#  \f$E(z; \cdots)\f$ is the [dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq),
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>      -#  \f$z\f$ is the redshift,
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{dV_C}{V_H}\f$, that is, **normalized** to the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3).<br>
    !>  To obtain the **full-sky (normalized) Comoving Volume Element**, multiply the output of the procedures under this generic interface by \f$4\pi\f$, that is, `4 * acos(-1.)`.<br>
    !>
    !>  \param[out] volComDiffNormed    :   The output scalar or array of the same rank as other array-like arguments,
    !>                                      of the same type and kind as the input argument `zplus1` containing the cosmological Comoving Volume Element
    !>                                      per [Steradian](https://en.wikipedia.org/wiki/Steradian) at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>  \param[in]  disComTransNormedSq :   The input scalar or array of the same rank as other array-like arguments, of type `real` of kind \RKALL representing the
    !>                                      [the square of the Dimensionless (Normalized) Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed) at the desired redshift for which the **Comoving Volume Element** must be computed.<br>
    !>                                      This argument can be readily obtained by **taking the square** of the output of [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed),
    !>                                      or if the cosmology is the Concordance model, one can use the approximation provided by [getDisComTransNormedWU10](@ref pm_cosmology::getDisComTransNormedWU10).<br>
    !>                                      See the examples below for usage.<br>
    !>  \param[in]  hubbleParamNormed   :   The input scalar or array of the same rank as other array-like arguments, of the same type and kind as `disComTransNormedSq` representing the
    !>                                      [Dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq) at the desired redshift for which the **Comoving Volume Element** must be computed.<br>
    !>                                      This argument can be readily obtained by **taking the square-root** of the output of [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq).<br>
    !>                                      See the examples below for usage.<br>
    !>
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: setVolComDiffNormed
    !>
    !>      call setVolComDiffNormed(volComDiffNormed, disComTransNormedSq, hubbleParamNormed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input arguments must be computed for the **same redshift and cosmological parameters**.<br>
    !>  The equivalence and consistencies of the input arguments are **not** verified within the algorithm
    !>  because such validations require extra information that is not provided as input.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The reason for requesting the rather awkward input arguments above instead of the redshift and other primitive cosmological parameters is to
    !>  ensure the computation efficiency of the algorithms at the highest level by taking out the potentially unnecessary computations out of the algorithms.<br>
    !>  See [getVolComDiffNormed](@ref pm_cosmology::getVolComDiffNormed) for an equivalent, more flexible, but potentially less performant interface.<br>
    !>
    !>  \devnote
    !>  There is no performance benefit in passing the inverse of `hubbleParamNormed` instead of the what is passed in the current interface.<br>
    !>  Do not attempt to change the interface.<br>
    !>
    !>  \see
    !>  [getVolComNormed](@ref pm_cosmology::getVolComNormed)<br>
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/setVolComDiffNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/setVolComDiffNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/setVolComDiffNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/setVolComDiffNormed/setVolComDiffNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface setVolComDiffNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module subroutine setVolComDiffNormed_D0_RK5(volComDiffNormed, disComTransNormedSq, hubbleParamNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolComDiffNormed_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)                   :: volComDiffNormed
        real(RKG)   , intent(in)                    :: disComTransNormedSq, hubbleParamNormed
    end subroutine
#endif

#if RK4_ENABLED
    pure elemental module subroutine setVolComDiffNormed_D0_RK4(volComDiffNormed, disComTransNormedSq, hubbleParamNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolComDiffNormed_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)                   :: volComDiffNormed
        real(RKG)   , intent(in)                    :: disComTransNormedSq, hubbleParamNormed
    end subroutine
#endif

#if RK3_ENABLED
    pure elemental module subroutine setVolComDiffNormed_D0_RK3(volComDiffNormed, disComTransNormedSq, hubbleParamNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolComDiffNormed_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)                   :: volComDiffNormed
        real(RKG)   , intent(in)                    :: disComTransNormedSq, hubbleParamNormed
    end subroutine
#endif

#if RK2_ENABLED
    pure elemental module subroutine setVolComDiffNormed_D0_RK2(volComDiffNormed, disComTransNormedSq, hubbleParamNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolComDiffNormed_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)                   :: volComDiffNormed
        real(RKG)   , intent(in)                    :: disComTransNormedSq, hubbleParamNormed
    end subroutine
#endif

#if RK1_ENABLED
    pure elemental module subroutine setVolComDiffNormed_D0_RK1(volComDiffNormed, disComTransNormedSq, hubbleParamNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVolComDiffNormed_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)                   :: volComDiffNormed
        real(RKG)   , intent(in)                    :: disComTransNormedSq, hubbleParamNormed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **full-sky** (\f$4\pi\f$ [Steradian](https://en.wikipedia.org/wiki/Steradian)) cosmological <b>Comoving Volume</b>
    !>  *normalized* to [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3), given the user-specified cosmological parameters.<br>
    !>
    !>  \details
    !>  The comoving volume \f$V_C\f$ is the volume measure in which number densities of non-evolving objects locked into Hubble Flow are constant with redshift.<br>
    !>  It is the proper volume times three factors of the relative scale factor now to then, or \f$(1 + z)^3\f$.<br>
    !>  The derivative of comoving distance with redshift is \f$\frac{dD_C}{dz} = \frac{1}{E(z)}\f$ where \f$E(z)\f$ is the [Dimensionless Hubble Parameter](@ref pm_cosmology::getHubbleParamNormedSq).<br>
    !>  Assuming \f$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)\f$ represent the normalized densities of Dark Matter, Dark Energy, Radiation Energy, and Curvature
    !>  in a Universe with negligible neutrino mass such that \f$\Omega_M + \Omega_\Lambda + \Omega_R + \Omega_K = 1\f$, the **all-sky** Comoving Volume \f$V_C\f$
    !>  out to a given cosmological redshift \f$z\f$ is defined as (e.g., Peebles, 1993),
    !>
    !>  \f{equation}{
    !>      \large
    !>      V_C(z; \Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K)
    !>      \begin{cases}
    !>          \bigg(\frac{4\pi V_H}{2\Omega_K}\bigg) ~ \bigg[ \frac{D_M}{D_H} \sqrt{1 + \Omega_K \bigg(\frac{D_M}{D_H}\bigg)^2} - \frac{1}{\sqrt{|\Omega_K|}}\mathrm{arcsinh}\bigg(\sqrt{|\Omega_K|}\frac{D_M}{D_H}\bigg) \bigg] &, \Omega_K > 0 ~, \\
    !>          \bigg(\frac{4\pi V_H}{2\Omega_K}\bigg) ~ \bigg[ \frac{D_M}{D_H} \sqrt{1 + \Omega_K \bigg(\frac{D_M}{D_H}\bigg)^2} - \frac{1}{\sqrt{|\Omega_K|}}\arcsin\bigg(\sqrt{|\Omega_K|}\frac{D_M}{D_H}\bigg) \bigg] &, \Omega_K < 0 ~, \\
    !>          \frac{4\pi}{3} D_M^3 &, \Omega_K = 0 ~,
    !>      \end{cases}
    !>  \f}
    !>
    !>  where
    !>      -#  \f$V_C(z; \cdots)\f$ is the cosmological Comoving Volume out to redshift \f$z\f$,
    !>      -#  \f$D_M(z; \cdots)\f$ is the cosmological [Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed) as a function of redshift,
    !>      -#  \f$\Omega_K\f$ is the [Curvature Density](@ref pm_cosmology::OMEGA_K),
    !>      -#  \f$D_H = \frac{C}{H_0}\f$ is the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC),
    !>      -#  \f$V_H = D_H^3\f$ is the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3),
    !>      -#  \f$H_0\f$ is the [Hubble Constant](@ref pm_cosmology::HUBBLE_CONST).
    !>      -#  \f$C\f$ is the [speed of light](@ref pm_cosmology::LIGHT_SPEED),
    !>      -#  \f$z\f$ is the redshift,
    !>
    !>  **The value returned by the procedures under this generic interface is** \f$\frac{V_C}{V_H}\f$, that is, **normalized** to the [Hubble Volume](@ref pm_cosmology::HUBBLE_VOLUME_MPC3).<br>
    !>
    !>  \param[in]  disComTransNormed   :   The input scalar or array of the same rank as other array-like arguments, of type `real` of kind \RKALL representing the
    !>                                      [Dimensionless (Normalized) Transverse Comoving Distance](@ref pm_cosmology::getDisComTransNormed) at the desired redshift for which the full-sky **Comoving Volume** must be computed.<br>
    !>                                      This argument can be readily obtained by **taking the square** of the output of [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed).<br>
    !>                                      See the examples below for usage.<br>
    !>  \param[in]  omegaK              :   The input scalar or array of the same rank as other array-like arguments,
    !>                                      of the same type and kind as the input argument `disComTransNormed` representing the normalized curvature density of the universe.<br>
    !>                                      (**optional**, default = [OMEGA_K](@ref pm_cosmology::OMEGA_K). It must be present <b>if and only if</b> `sqrtAbsOmegaK` is also present.)
    !>  \param[in]  sqrtAbsOmegaK       :   The input scalar or array of the same rank as other array-like arguments,
    !>                                      of the same type and kind as the input argument `disComTransNormed` representing the square root of the normalized curvature density of the universe `sqrt(omegaK)`.<br>
    !>                                      This argument is requested as an input to avoid the redundant costly `sqrt(OmegaK)` operation within the algorithm
    !>                                      when the procedure is to be called repeatedly with the same `omegaK`.<br>
    !>                                      (**optional**, it must be present <b>if and only if</b> `omegaK` is also present.)
    !>
    !>  \return
    !>  `volComNormed`                  :   The output scalar or array of the same rank as other array-like arguments,
    !>                                      of the same type and kind as the input argument `disComTransNormed` containing the cosmological full-sky Comoving Volume
    !>                                      at the desired redshift **normalized** to the [Hubble Distance](@ref pm_cosmology::HUBBLE_DISTANCE_MPC).
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_cosmology, only: getVolComNormed
    !>
    !>      volComNormed = getVolComNormed(disComTransNormed)
    !>      volComNormed = getVolComNormed(disComTransNormed, omegaK, sqrtAbsOmegaK)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input arguments must be computed for the **same redshift and cosmological parameters**.<br>
    !>  The equivalence and consistencies of the input arguments are **not** verified within the algorithm because such validations require extra information that is not provided as input.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  The reason for requesting the rather awkward input arguments above instead of the redshift and other primitive cosmological parameters is to
    !>  ensure the computation efficiency of the algorithms at the highest level by taking out the potentially unnecessary computations out of the algorithms.<br>
    !>
    !>  \note
    !>  Dropping all optional arguments corresponds to a flat Universe.<br>
    !>
    !>  \see
    !>  [getVolComDiffNormed](@ref pm_cosmology::getVolComDiffNormed)<br>
    !>  [getHubbleParamNormedSq](@ref pm_cosmology::getHubbleParamNormedSq)<br>
    !>  [getDisComTransNormed](@ref pm_cosmology::getDisComTransNormed)<br>
    !>  [getDisLumNormed](@ref pm_cosmology::getDisLumNormed)<br>
    !>  [getDisAngNormed](@ref pm_cosmology::getDisAngNormed)<br>
    !>  [getDisComNormed](@ref pm_cosmology::getDisComNormed)<br>
    !>  [LOG_HUBBLE_CONST](@ref pm_cosmology::LOG_HUBBLE_CONST)<br>
    !>  [HUBBLE_DISTANCE_MPC](@ref pm_cosmology::HUBBLE_DISTANCE_MPC)<br>
    !>  [HUBBLE_CONST](@ref pm_cosmology::HUBBLE_CONST)<br>
    !>  [OMEGA_M](@ref pm_cosmology::OMEGA_M)<br>
    !>  [OMEGA_L](@ref pm_cosmology::OMEGA_L)<br>
    !>  [OMEGA_R](@ref pm_cosmology::OMEGA_R)<br>
    !>  [OMEGA_K](@ref pm_cosmology::OMEGA_K)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_cosmology/getVolComNormed/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_cosmology/getVolComNormed/main.out.F90
    !>  \postproc
    !>  \include{lineno} example/pm_cosmology/getVolComNormed/main.py
    !>  \vis
    !>  \image html pm_cosmology/getVolComNormed/getVolComNormed.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_cosmology](@ref test_pm_cosmology)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    interface getVolComNormed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getVolComNormedD_D0_RK5(disComTransNormed) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedD_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: disComTransNormed
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getVolComNormedD_D0_RK4(disComTransNormed) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedD_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: disComTransNormed
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getVolComNormedD_D0_RK3(disComTransNormed) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedD_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: disComTransNormed
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getVolComNormedD_D0_RK2(disComTransNormed) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedD_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: disComTransNormed
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getVolComNormedD_D0_RK1(disComTransNormed) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedD_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: disComTransNormed
        real(RKG)                                   :: volComNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getVolComNormedDOS_D0_RK5(disComTransNormed, omegaK, sqrtAbsOmegaK) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedDOS_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)                    :: disComTransNormed, omegaK, sqrtAbsOmegaK
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getVolComNormedDOS_D0_RK4(disComTransNormed, omegaK, sqrtAbsOmegaK) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedDOS_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)                    :: disComTransNormed, omegaK, sqrtAbsOmegaK
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getVolComNormedDOS_D0_RK3(disComTransNormed, omegaK, sqrtAbsOmegaK) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedDOS_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)                    :: disComTransNormed, omegaK, sqrtAbsOmegaK
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getVolComNormedDOS_D0_RK2(disComTransNormed, omegaK, sqrtAbsOmegaK) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedDOS_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)                    :: disComTransNormed, omegaK, sqrtAbsOmegaK
        real(RKG)                                   :: volComNormed
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getVolComNormedDOS_D0_RK1(disComTransNormed, omegaK, sqrtAbsOmegaK) result(volComNormed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVolComNormedDOS_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)                    :: disComTransNormed, omegaK, sqrtAbsOmegaK
        real(RKG)                                   :: volComNormed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

#if LEGACY_ENABLED
! LCOV_EXCL_START
!>  \cond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Return the approximation to the cosmological luminosity distance.
    !>
    !>  \param[in]   zplus1 : The redshift plus 1.
    !>
    !>  \return
    !>  `ldiswickram` : The approximate cosmological luminosity distance.
    !>
    !>  The approximation is accurate with a relative error of 0.001 for any redshift above 0.1.
    !>
    !>  \remark
    !>  This function is the same as [getLogLumDisWicMpc()](@ref getLogLumDisWicMpc) function,
    !>  except for the fact that it return the natural value, as opposed to the natural logarithm.
    !>  It is kept only for legacy reasons and should not be used in new code.
    !>
    !>  \author
    !>  \AmirShahmoradi, Sunday 2:31 PM, January 6, 2013, Institute for Fusion Studies, The University of Texas at Austin.<br>
    pure elemental function ldiswickram(zplus1)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ldiswickram
#endif
        use pm_kind, only: RK
        real(RK), intent(in)    :: zplus1
        real(RK)                :: ldiswickram,alpha,alpha0,x,x0,psi,psi0
        real(RK), parameter     :: OMEGA_M_RK = real(OMEGA_M, RK)
        real(RK), parameter     :: OMEGA_L_RK = real(OMEGA_L, RK)
        alpha = 1._RK + 2 * OMEGA_L_RK / (OMEGA_M_RK*zplus1**3)
        alpha0 = 1._RK + 2 * OMEGA_L_RK / OMEGA_M_RK
        x0 = log(alpha0+sqrt(alpha0*alpha0 - 1._RK))
        x = log(alpha + sqrt(alpha * alpha - 1._RK))
        psi=x**0.333333333333333_RK*(1.5874010519682_RK-6.2992105236833e-3_RK*x*x+7.5375168659459e-5_RK*x**4)
        psi0=x0**0.333333333333333_RK*(1.5874010519682_RK-6.2992105236833e-3_RK*x0*x0+7.5375168659459e-5_RK*x0**4)
        ldiswickram = real(HUBBLE_DISTANCE_MPC, RK) * zplus1 * (psi0 - psi) / (OMEGA_L_RK**0.1666666666666667_RK * OMEGA_M_RK**0.333333333333333_RK)
    end function ldiswickram

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Return the cosmological lookback time in GYrs at the given redshift for the assumed cosmological parameters.
    !>
    !>  \param[in]   zplus1              : The redshift plus 1.
    !>  \param[in]   maxRelativeError    : The maximum tolerance for error in the numerical integration (**optional**, default = 1.e-6).
    !>  \param[in]   nRefinement         : The number of refinements in the Romberg numerical integration (**optional**, default = 5).
    !>  \return
    !>  `lookBackTime` : The cosmological lookback time in GYrs at the given redshift.
    !>
    !>  \remark
    !>  The rest of the required parameters are taken from the module:
    !>      - `HUBBLE_TIME_GYR` : The Hubble time in units of Gyrs (Giga Years). It should be a number close to 13.8 Gyrs (Liddle, 2003, Page 57).
    !>
    !>  \remark
    !>  The integrations are performed using the Romberg integration method.
    !>
    !>  \author
    !>  \AmirShahmoradi, Sunday 2:31 PM, January 6, 2013, Institute for Fusion Studies, The University of Texas at Austin.<br>
    impure function getLookBackTime(zplus1, maxRelativeError, nRefinement) result(lookBackTime)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLookBackTime
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use pm_quadRomb, only: getQuadRomb
        use pm_kind, only: RK
        real(RK)    , intent(in)            :: zplus1
        real(RK)    , intent(in), optional  :: maxRelativeError
        integer(IK) , intent(in), optional  :: nRefinement
        real(RK)    , parameter             :: ZPLUS1_MIN = 1._RK
        real(RK)                            :: lookBackTime
        real(RK)                            :: maxRelativeErrorDefault
        integer(IK)                         :: nRefinementDefault

        maxRelativeErrorDefault = 1.e-6_RK; if (present(maxRelativeError)) maxRelativeErrorDefault = maxRelativeError
        nRefinementDefault = 7_IK; if (present(nRefinement)) nRefinementDefault = nRefinement

        lookBackTime = getQuadRomb  ( getFunc   = getLookBackTimeDensity    &
                                    , lb        = ZPLUS1_MIN                &
                                    , ub        = zplus1                    &
                                    , tol       = maxRelativeErrorDefault   &
                                    , nref      = nRefinementDefault        &
                                    )
        lookBackTime = real(HUBBLE_TIME_GYR, RK) * lookBackTime

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Return the **differential** (w.r.t. `z`) cosmological lookback time in GYrs at the given redshift for the assumed cosmological parameters.
    !>
    !>  \param[in]   zplus1              : The redshift plus 1.
    !>  \return
    !>  `lookBackTimeDnesity` : The cosmological lookback time in GYrs at the given redshift.
    !>
    !>  \remark
    !>  The rest of the required parameters are taken from the module:
    !>  - `OMEGA_L` : the dark energy density,
    !>  - `OMEGA_M` : the dark matter density.
    !>
    !>  \remark
    !>  This function could have become an internal function within [getLookBackTime()](@ref getLookBackTime).
    !>  However, the library yields segmentation fault error when compiled and run on the Windows Subsystem
    !>  for Linux Ubuntu with GFortran. As such, it is implemented as an independent function.
    !>
    !>  \author
    !>  \AmirShahmoradi, Sunday 2:31 PM, January 6, 2013, Institute for Fusion Studies, The University of Texas at Austin.<br>
    pure function getLookBackTimeDensity(zplus1) result(lookBackTimeDnesity)
        use pm_kind, only: RK
        real(RK), intent(in)    :: zplus1
        real(RK)                :: lookBackTimeDnesity
        real(RK), parameter     :: OMEGA_M_RK = real(OMEGA_M, RK)
        real(RK), parameter     :: OMEGA_L_RK = real(OMEGA_L, RK)
        lookBackTimeDnesity = 1._RK / ( zplus1 * sqrt(OMEGA_M_RK * zplus1**3 + OMEGA_L_RK) )
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Return the derivative of the age of the Universe, w.r.t. redshift for a given input redshift + 1.
    !>
    !>  \param[in] zplus1 : The redshift plus 1.
    !>
    !>  \return
    !>  `universeAgeDerivative` : The derivative of the age of the Universe, w.r.t. redshift for a given input `zplus1`.
    !>
    !>  \remark
    !>  The rest of the required parameters are taken from the module:
    !>  - `INV_HUBBLE_CONST`
    !>  - `OMEGA_L`
    !>  - `OMEGA_M`
    !>
    !>  \remark
    !>  The integrations are performed using the Romberg integration method.
    !>
    !>  \author
    !>  \AmirShahmoradi, Sunday 2:31 PM, January 6, 2013, Institute for Fusion Studies, The University of Texas at Austin.<br>
    pure elemental function getUniverseAgeDerivative(zplus1) result(universeAgeDerivative)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniverseAgeDerivative
#endif
        use pm_kind, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: universeAgeDerivative
        real(RK), parameter     :: OMEGA_M_RK = real(OMEGA_M, RK)
        real(RK), parameter     :: OMEGA_L_RK = real(OMEGA_L, RK)
        real(RK), parameter     :: INV_HUBBLE_CONST_RK = real(INV_HUBBLE_CONST, RK)
        universeAgeDerivative = INV_HUBBLE_CONST_RK / ( zplus1 * sqrt(OMEGA_M_RK * zplus1**3 + OMEGA_L_RK) )
    end function getUniverseAgeDerivative

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \endcond excluded
! LCOV_EXCL_STOP
#endif

end module pm_cosmology