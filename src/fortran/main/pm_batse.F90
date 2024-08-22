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
!>  This module contains procedures and generic interfaces for modeling data and
!>  detectors of the BATSE Gamma-Ray detectors onboard the NASA Compton Gamma-Ray Observatory.<br>
!>
!>  \details
!>  The Compton Gamma Ray Observatory (CGRO) was a space observatory detecting photons with energies from 20 keV to 30 GeV, in Earth orbit from 1991 to 2000.<br>
!>  The observatory featured four main telescopes in one spacecraft, covering X-rays and gamma rays, including various specialized sub-instruments and detectors.<br>
!>  Following 14 years of effort, the observatory was launched from Space Shuttle Atlantis during STS-37 on April 5, 1991, and operated until its deorbit on June 4, 2000.<br>
!>  It was deployed in low Earth orbit at 450 km (280 mi) to avoid the Van Allen radiation belt.<br>
!>  It was the heaviest astrophysical payload ever flown at that time at 16,300 kilograms (35,900 lb).<br>
!>
!>  Instruments
!>  ===========
!>
!>  CGRO carried a complement of four instruments that covered an unprecedented six orders of the electromagnetic spectrum, from 20 keV to 30 GeV.<br>
!>  Only the major instrument of interest to this module is discussed below.<br>
!>
!>  BATSE
!>  -----
!>
!>  \image html pm_batse.gif width=500
!>
!>  The Burst and Transient Source Experiment (BATSE) by the NASA Marshall Space Flight Center searched
!>  the sky for gamma-ray bursts (20 to >600 keV) and conducted full-sky surveys for long-lived sources.<br>
!>  It consisted of eight identical detector modules, one at each of the satellite corners.<br>
!>  Each module consisted of both a NaI(Tl) Large Area Detector (LAD) covering the 20 keV to ~2 MeV range, 50.48 cm in dia by 1.27 cm thick,
!>  and a 12.7 cm dia by 7.62 cm thick NaI Spectroscopy Detector, which extended the upper energy range to 8 MeV, all surrounded by a
!>  plastic scintillator in active anti-coincidence to veto the large background rates due to cosmic rays and trapped radiation.<br>
!>  Sudden increases in the LAD rates triggered a high-speed data storage mode, the details of the burst being read out to telemetry later.<br>
!>  Bursts were typically detected at rates of roughly one per day over the 9-year CGRO mission.<br>
!>  A strong burst could result in the observation of many thousands of gamma-rays within a time interval ranging from ~0.1 s up to about 100 s.<br>
!>
!>  Other instruments
!>  -----------------
!>
!>  See the references below.<br>
!>
!>  \see
!>  [pm_distBand](@ref pm_distBand)<br>
!>  Kaneko, 2005, Spectral studies of gamma-ray burst prompt emission.<br>
!>  [Shahmoradi and Nemiroff, 2015, MNRAS, 451:4645-4662.](https://www.cdslab.org/pubs/2015_Shahmoradi_I.pdf)<br>
!>
!>  \test
!>  [test_pm_batse](@ref test_pm_batse)<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Tuesday, April 30, 2019, 12:58 PM, SEIR, UTA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_batse

    use pm_kind, only: SK, IK, RKD, RKB

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_batse"

    real(RKB)       , parameter :: LOG_TEN = log(10._RKB)
    real(RKB)       , parameter :: MIN_LOG10PH53_4_LOGPBOLZERO  = 4.92_RKB
    real(RKB)       , parameter :: MAX_LOG10PH53_4_LOGPBOLZERO  = 6.318167895318538_RKB
    real(RKB)       , parameter :: MAX_LOGPH53_4_LOGPBOLZERO    = MAX_LOG10PH53_4_LOGPBOLZERO * LOG_TEN
    real(RKB)       , parameter :: MIN_LOGPH53_4_LOGPBOLZERO    = MIN_LOG10PH53_4_LOGPBOLZERO * LOG_TEN
    real(RKB)       , parameter :: DIF_LOGPH53_4_LOGPBOLZERO    = MAX_LOGPH53_4_LOGPBOLZERO - MIN_LOGPH53_4_LOGPBOLZERO

    ! Parameters for the effective peak flux of SGRBs (P64ms conversion to P1024):

!   ! The scale of the change in BATSE efficiency for different GRB durations.
!   real(RKB)       , parameter :: THRESH_ERFC_BASE             = +0.146314238936889_RKB
!
!   ! The scale of the change in BATSE efficiency for different GRB durations.
!   real(RKB)       , parameter :: THRESH_ERFC_AMP              = +0.570285374263156_RKB
!
!   ! Mean duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
!   real(RKB)       , parameter :: THRESH_ERFC_AVG              = -0.480811530417719_RKB
!
!   ! scale of the duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
!   real(RKB)       , parameter :: THRESH_ERFC_STD              = +1.292443569094922_RKB

    real(RKB)       , parameter :: LOGPF53_MINUS_LOGPBOL        = 11.328718657530706_RKB
    real(RKB)       , parameter :: LOG10PF53_MINUS_LOG10PBOL    = LOGPF53_MINUS_LOGPBOL / LOG_TEN

    !>  The scale of the change in BATSE efficiency for different SGRB durations.
    real(RKB)       , parameter :: THRESH_ERFC_BASE             = +0.146314238936889_RKB * LOG_TEN

    !>  The scale of the change in BATSE efficiency for different GRB durations.
    real(RKB)       , parameter :: THRESH_ERFC_AMP              = +0.282313526464596_RKB * LOG_TEN ! Serfc

    !>  Mean duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
    real(RKB)       , parameter :: THRESH_ERFC_AVG              = -0.483553339256463_RKB * LOG_TEN  ! meandur

    !>  Scale of the duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
    real(RKB)       , parameter :: THRESH_ERFC_STD              = 1.0514698984694800_RKB * LOG_TEN  ! scaledur

    !>  Inverse scale of the duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
    real(RKB)       , parameter :: THRESH_ERFC_STD_INV          = 1._RKB / THRESH_ERFC_STD       ! inverse scaledur

    ! The height of the ERFC function.
    real(RKB)       , parameter :: THRESH_ERFC_HEIGHT           = 2 * THRESH_ERFC_AMP

    !>  Correction that must be added to logPbol64ms to convert it to effective peak flux.
    !   Effective LogPbol limit above which trigger efficiency is 100%, for any Log(Epk) and Log(dur). It is equivalent to maximum Log(Pbol) at very long durations.
    real(RKB)       , parameter :: THRESH_LOGPBOL64_CORRECTION  = DIF_LOGPH53_4_LOGPBOLZERO - MIN_LOGPH53_4_LOGPBOLZERO + THRESH_ERFC_HEIGHT ! equivalent to lpb_correction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating objects that contain attributes of BATSE catalog GRBs.<br>
    !>
    !>  \test
    !>  [test_pm_batse](@ref test_pm_batse)
    !>
    !>  \final{getCorrectionLogEffPPF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type :: grb_type
        real(RKD) :: pbol    = -huge(0._RKD) ! < \public The scalar `real` of kind \RKD, containing a measure of the bolometric energy flux of a BATSE GRB event.
        real(RKD) :: epeak   = -huge(0._RKD) ! < \public The scalar `real` of kind \RKD, containing a measure of the spectral peak energy of a BATSE GRB event.
        real(RKD) :: sbol    = -huge(0._RKD) ! < \public The scalar `real` of kind \RKD, containing a measure of the bolometric energy fluence of a BATSE GRB event.
        real(RKD) :: t90     = -huge(0._RKD) ! < \public The scalar `real` of kind \RKD, containing a measure of the duration of a BATSE GRB event.
        real(RKD) :: pph1024 = -huge(0._RKD) ! < \public The scalar `real` of kind \RKD, containing a measure of the peak 1024ms photon flux of a BATSE GRB event in the BATSE nominal detection energy window \f$[50, 300]\f$.
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the correction required for converting an input **natural-log** peak photon flux in 64ms timescale in
    !>  the BATSE detection energy range to an effective triggering peak photon flux in the same energy range.<br>
    !>
    !>  \details
    !>  See the *Eqn A4 of Shahmoradi and Nemiroff 2015, MNRAS, Short versus long gamma-ray bursts* for more details.<br>
    !>  The observed T90 duration of the event is necessary as input to account for varying duration of short GRBs.<br>
    !>
    !>  \param[in]  logT90  :   The input scalar, or array of arbitrary rank, of type `real` of kind \RKALL,
    !>                          containing the natural logarithm of the T90 duration-definition measure of a GRB in units of seconds.<br>
    !>
    !>  \return
    !>  `correction`        :   The output scalar or array of the same rank as the input argument of same type and kind as the input argument,
    !>                          containing the correction that must be subtracted from the natural logarithm of a given 64ms peak flux
    !>                          to convert it to an effective 1024ms triggering peak flux.
    !>
    !>  \interface{getCorrectionLogEffPPF}
    !>  \code{.F90}
    !>
    !>      use pm_batse, only: getCorrectionLogEffPPF
    !>
    !>      correction = getCorrectionLogEffPPF(logT90)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The following relationship holds: `logPPF64 - getCorrectionLogEffPPF(logT90) = logPPF1024 - THRESH_ERFC_BASE`.<br>
    !>
    !>  \see
    !>  [getLogPbol](@ref pm_batse::getLogPbol)<br>
    !>  [getLogPF53](@ref pm_batse::getLogPF53)<br>
    !>  [getLog10PF53](@ref pm_batse::getLog10PF53)<br>
    !>  [getLogEffPPF](@ref pm_batse::getLogEffPPF)<br>
    !>  [getCorrectionLogEffPPF](@ref pm_batse::getCorrectionLogEffPPF)<br>
    !>
    !>  \test
    !>  [test_pm_batse](@ref test_pm_batse)
    !>
    !>  \final{getCorrectionLogEffPPF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getCorrectionLogEffPPF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getCorrectionLogEffPPF_RK5(logT90) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorrectionLogEffPPF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: logT90
        real(RKG)                   :: correction
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getCorrectionLogEffPPF_RK4(logT90) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorrectionLogEffPPF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: logT90
        real(RKG)                   :: correction
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getCorrectionLogEffPPF_RK3(logT90) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorrectionLogEffPPF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: logT90
        real(RKG)                   :: correction
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getCorrectionLogEffPPF_RK2(logT90) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorrectionLogEffPPF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: logT90
        real(RKG)                   :: correction
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getCorrectionLogEffPPF_RK1(logT90) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorrectionLogEffPPF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: logT90
        real(RKG)                   :: correction
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getCorrectionLogEffPPF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of an input natural logarithm of peak photon flux in 64ms timescale to
    !>  an effective triggering peak photon flux both in the BATSE detection energy range \f$[50, 300]\kev\f$.<br>
    !>
    !>  \details
    !>  See the *Eqn A4 of Shahmoradi and Nemiroff 2015, MNRAS, Short versus long gamma-ray bursts* for more details.<br>
    !>  The observed T90 duration of the event is necessary as input to account for varying duration of short GRBs and their effects on the triggering mechanism.<br>
    !>
    !>  \param[in]  logPPF64    :   The input scalar, or array of the same rank as other input array-valued arguments, of type `real` of kind \RKALL,
    !>                              containing the natural logarithm of the 64ms Peak Photon Flux of an event in units of \f$[\ms{photons} / \ms{sec}]\f$.<br>
    !>  \param[in]  logT90      :   The input scalar, or array of the same rank as other input array-valued arguments, of the same type and kind as `logPPF64`,
    !>                              containing the natural logarithm of the T90 duration-definition measure of a GRB in units of seconds.<br>
    !>
    !>  \return
    !>  `logEffPPF`             :   The output scalar or array of the same rank as the input arguments of same type and kind as the input arguments,
    !>                              containing the natural logarithm of the effective 1024ms triggering peak photon flux in the BATSE nominal
    !>                              detection energy range, such that the effects of multiple triggering timescale of BATSE detectors
    !>                              are compensated for, as if all BATSE GRBs with different observed durations were detected
    !>                              by the same duration-agnostic gamma-ray detector and triggering mechanism.
    !>
    !>  \interface{getLogEffPPF}
    !>  \code{.F90}
    !>
    !>      use pm_batse, only: getLogEffPPF
    !>
    !>      logEffPPF = getLogEffPPF(logPPF64, logT90)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The following relationship holds: `logPPF64 - getCorrectionLogEffPPF(logT90) = logPPF1024 - THRESH_ERFC_BASE`.<br>
    !>
    !>  \see
    !>  [getLogPbol](@ref pm_batse::getLogPbol)<br>
    !>  [getLogPF53](@ref pm_batse::getLogPF53)<br>
    !>  [getLog10PF53](@ref pm_batse::getLog10PF53)<br>
    !>  [getLogEffPPF](@ref pm_batse::getLogEffPPF)<br>
    !>  [getCorrectionLogEffPPF](@ref pm_batse::getCorrectionLogEffPPF)<br>
    !>
    !>  \test
    !>  [test_pm_batse](@ref test_pm_batse)
    !>
    !>  \final{getLogEffPPF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getLogEffPPF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getLogEffPPF_RK5(logPPF64, logT90) result(logEffPPF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEffPPF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: logPPF64, logT90
        real(RKG)                   :: logEffPPF
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getLogEffPPF_RK4(logPPF64, logT90) result(logEffPPF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEffPPF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: logPPF64, logT90
        real(RKG)                   :: logEffPPF
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getLogEffPPF_RK3(logPPF64, logT90) result(logEffPPF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEffPPF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: logPPF64, logT90
        real(RKG)                   :: logEffPPF
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getLogEffPPF_RK2(logPPF64, logT90) result(logEffPPF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEffPPF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: logPPF64, logT90
        real(RKG)                   :: logEffPPF
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getLogEffPPF_RK1(logPPF64, logT90) result(logEffPPF)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEffPPF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: logPPF64, logT90
        real(RKG)                   :: logEffPPF
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLogEffPPF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of an input natural logarithm of photon flux/fluence over a given time in the BATSE detection energy
    !>  range \f$[50, 300]\kev\f$ to a bolometric (\f$[0.0001, 20000]\kev\f$) energy flux/fluence in units of \f$\ergs\f$ over the same time.<br>
    !>
    !>  \details
    !>  See the *Eqn A4 of Shahmoradi and Nemiroff 2015, MNRAS, Short versus long gamma-ray bursts* for more details.<br>
    !>
    !>  \param[in]  logEpk      :   The input scalar, or array of the same rank as other input array-valued arguments, of the same type and kind as `logPPF64`,
    !>                              containing the natural logarithm of the spectral peak energy of the GRB in units of \f$\ergs\f$.<br>
    !>  \param[in]  logPF53     :   The input scalar, or array of the same rank as other input array-valued arguments, of type `real` of kind \RKALL,
    !>                              containing the natural logarithm of the Photon Flux of an event in units of \f$[\ms{photons} / \ms{sec}]\f$.<br>
    !>
    !>  \return
    !>  `logPbol`               :   The output scalar or array of the same rank as the input arguments of same type and kind as the input arguments,
    !>                              containing the natural logarithm of the bolometric (\f$[0.0001, 20000]\kev\f$) energy flux of the GRB event
    !>                              in units of \f$\ergs\f$ over the same time interval.<br>
    !>
    !>  \interface{getLogPbol}
    !>  \code{.F90}
    !>
    !>      use pm_batse, only: getLogPbol
    !>
    !>      logPbol = getLogPbol(logEpk, logPF53)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogPbol](@ref pm_batse::getLogPbol)<br>
    !>  [getLogPF53](@ref pm_batse::getLogPF53)<br>
    !>  [getLog10PF53](@ref pm_batse::getLog10PF53)<br>
    !>  [getLogEffPPF](@ref pm_batse::getLogEffPPF)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [getCorrectionLogEffPPF](@ref pm_batse::getCorrectionLogEffPPF)<br>
    !>
    !>  \example{getLogPbol}
    !>  \include{lineno} example/pm_batse/getLogPbol/main.F90
    !>  \compilef{getLogPbol}
    !>  \output{getLogPbol}
    !>  \include{lineno} example/pm_batse/getLogPbol/main.out.F90
    !>  \postproc{getLogPbol}
    !>  \include{lineno} example/pm_batse/getLogPbol/main.py
    !>  \vis{getLogPbol}
    !>  \image html pm_batse/getLogPbol/getLogPbol.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_batse](@ref test_pm_batse)
    !>
    !>  \final{getLogPbol}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday June 27, 2012, 7:15 PM, Institute for Fusion Studies, The University of Texas at Austin.<br>
    interface getLogPbol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getLogPbol_RK5(logEpk, logPF53) result(logPbol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPbol_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: logEpk, logPF53
        real(RKG)                   :: logPbol
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getLogPbol_RK4(logEpk, logPF53) result(logPbol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPbol_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: logEpk, logPF53
        real(RKG)                   :: logPbol
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getLogPbol_RK3(logEpk, logPF53) result(logPbol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPbol_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: logEpk, logPF53
        real(RKG)                   :: logPbol
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getLogPbol_RK2(logEpk, logPF53) result(logPbol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPbol_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: logEpk, logPF53
        real(RKG)                   :: logPbol
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getLogPbol_RK1(logEpk, logPF53) result(logPbol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPbol_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: logEpk, logPF53
        real(RKG)                   :: logPbol
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLogPbol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of an input natural logarithm of a bolometric (\f$[0.0001, 20000]\kev\f$) energy flux/fluence
    !>  in units of \f$\ergs\f$ over a given time to the photon flux/fluence over the same time in the BATSE detection energy range \f$[50, 300]\kev\f$.<br>
    !>
    !>  \details
    !>  See the *Eqn A4 of Shahmoradi and Nemiroff 2015, MNRAS, Short versus long gamma-ray bursts* for more details.<br>
    !>
    !>  \param[in]  logEpk      :   The input scalar, or array of the same rank as other input array-valued arguments, of the same type and kind as `logPPF64`,
    !>                              containing the natural logarithm of the spectral peak energy of the GRB in units of \f$\ergs\f$.<br>
    !>  \param[in]  logPbol     :   The input scalar or array of the same rank as the input arguments of same type and kind as `logPPF64`,
    !>                              containing the natural logarithm of the bolometric (\f$[0.0001, 20000]\kev\f$) energy flux of the GRB event
    !>                              in units of \f$\ergs\f$ over the same time interval.<br>
    !>
    !>  \return
    !>  `logPF53`               :   The output scalar, or array of the same rank as other input array-valued arguments, of type `real` of kind \RKALL,
    !>                              containing the natural logarithm of the Photon Flux of an event in units of \f$[\ms{photons} / \ms{sec}]\f$.<br>
    !>
    !>  \interface{getLogPF53}
    !>  \code{.F90}
    !>
    !>      use pm_batse, only: getLogPF53
    !>
    !>      logPF53 = getLogPF53(logEpk, logPbol)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogPbol](@ref pm_batse::getLogPbol)<br>
    !>  [getLogPF53](@ref pm_batse::getLogPF53)<br>
    !>  [getLog10PF53](@ref pm_batse::getLog10PF53)<br>
    !>  [getLogEffPPF](@ref pm_batse::getLogEffPPF)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [getCorrectionLogEffPPF](@ref pm_batse::getCorrectionLogEffPPF)<br>
    !>
    !>  \example{getLogPF53}
    !>  \include{lineno} example/pm_batse/getLogPF53/main.F90
    !>  \compilef{getLogPF53}
    !>  \output{getLogPF53}
    !>  \include{lineno} example/pm_batse/getLogPF53/main.out.F90
    !>  \postproc{getLogPF53}
    !>  \include{lineno} example/pm_batse/getLogPF53/main.py
    !>  \vis{getLogPF53}
    !>  \image html pm_batse/getLogPF53/getLogPF53.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_batse](@ref test_pm_batse)
    !>
    !>  \final{getLogPF53}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday June 27, 2012, 7:15 PM, Institute for Fusion Studies, The University of Texas at Austin.<br>
    interface getLogPF53

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getLogPF53_RK5(logEpk, logPbol) result(logPF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPF53_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: logEpk, logPbol
        real(RKG)                   :: logPF53
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getLogPF53_RK4(logEpk, logPbol) result(logPF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPF53_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: logEpk, logPbol
        real(RKG)                   :: logPF53
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getLogPF53_RK3(logEpk, logPbol) result(logPF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPF53_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: logEpk, logPbol
        real(RKG)                   :: logPF53
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getLogPF53_RK2(logEpk, logPbol) result(logPF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPF53_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: logEpk, logPbol
        real(RKG)                   :: logPF53
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getLogPF53_RK1(logEpk, logPbol) result(logPF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPF53_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: logEpk, logPbol
        real(RKG)                   :: logPF53
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLogPF53

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  This generic interface is identical to the generic interface [getLogPF53](@ref pm_batse::getLogPF53)
    !>  with the only difference that all input and output arguments must be in \f$\log_{10}(\cdot)\f$ instead of natural logarithm.
    !>
    !>  \details
    !>  See the documentation of [getLogPF53](@ref pm_batse::getLogPF53) for more information.<br>
    !>
    !>  \param[in]  log10Epk    :   The input scalar, or array of the same rank as other input array-valued arguments, of the same type and kind as `logPPF64`,
    !>                              containing the log10 of the spectral peak energy of the GRB in units of \f$\ergs\f$.<br>
    !>  \param[in]  log10Pbol   :   The input scalar or array of the same rank as the input arguments of same type and kind as `logPPF64`,
    !>                              containing the log10 of the bolometric (\f$[0.0001, 20000]\kev\f$) energy flux of the GRB event
    !>                              in units of \f$\ergs\f$ over the same time interval.<br>
    !>
    !>  \return
    !>  `log10PF53`             :   The output scalar, or array of the same rank as other input array-valued arguments, of type `real` of kind \RKALL,
    !>                              containing the log10 of the Photon Flux of an event in units of \f$[\ms{photons} / \ms{sec}]\f$.<br>
    !>
    !>  \interface{getLog10PF53}
    !>  \code{.F90}
    !>
    !>      use pm_batse, only: getLog10PF53
    !>
    !>      log10PF53 = getLog10PF53(log10Epk, log10Pbol)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLogPbol](@ref pm_batse::getLogPbol)<br>
    !>  [getLogPF53](@ref pm_batse::getLogPF53)<br>
    !>  [getLog10PF53](@ref pm_batse::getLog10PF53)<br>
    !>  [getLogEffPPF](@ref pm_batse::getLogEffPPF)<br>
    !>  [setBandEnergy](@ref pm_distBand::setBandEnergy)<br>
    !>  [setBandPhoton](@ref pm_distBand::setBandPhoton)<br>
    !>  [getCorrectionLogEffPPF](@ref pm_batse::getCorrectionLogEffPPF)<br>
    !>
    !>  \test
    !>  [test_pm_batse](@ref test_pm_batse)
    !>
    !>  \final{getLog10PF53}
    !>
    !>  \author
    !>  \AmirShahmoradi, Wednesday June 27, 2012, 7:15 PM, Institute for Fusion Studies, The University of Texas at Austin.<br>
    interface getLog10PF53

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getLog10PF53_RK5(log10Epk, log10Pbol) result(log10PF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog10PF53_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: log10Epk, log10Pbol
        real(RKG)                   :: log10PF53
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getLog10PF53_RK4(log10Epk, log10Pbol) result(log10PF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog10PF53_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: log10Epk, log10Pbol
        real(RKG)                   :: log10PF53
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getLog10PF53_RK3(log10Epk, log10Pbol) result(log10PF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog10PF53_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: log10Epk, log10Pbol
        real(RKG)                   :: log10PF53
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getLog10PF53_RK2(log10Epk, log10Pbol) result(log10PF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog10PF53_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: log10Epk, log10Pbol
        real(RKG)                   :: log10PF53
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getLog10PF53_RK1(log10Epk, log10Pbol) result(log10PF53)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog10PF53_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: log10Epk, log10Pbol
        real(RKG)                   :: log10PF53
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getLog10PF53

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !> Return all log of data in natural (Neper) base.
!    !>
!    !> \param[in]   inFilePath  :   The path to the input BATSE file.
!    !> \param[in]   outFilePath :   The path to the output BATSE file.
!    !> \param[in]   isLgrb      :   A logical flag indicating what type of input file is being processed.
!    ! integer(IK)     , parameter :: NLGRB = 1366_IK, NSGRB = 565_IK, NVAR = 4_IK
!
!    ! GRB attributes
!    type :: Event_type
!        real(RK) :: logPbol, logEpk, logSbol, logT90, logPF53
!    end type Event_type
!
!    type :: GRB_type
!        integer(IK) :: count
!        type(Event_type), allocatable :: Event(:)
!        !type(Event_type) :: Event(NLGRB)
!    end type GRB_type
!
!    !>  \cond excluded
!#if CAF_ENABLED
!    type(GRB_type) :: GRB[*]
!#else
!    type(GRB_type) :: GRB
!#endif
!    !>  \endcond excluded
!
!    integer(IK), allocatable :: Trigger(:)
!    !integer(IK) :: Trigger(NLGRB)
!    !integer(IK) :: TriggerSGRB(NSGRB)
!
!    subroutine readDataGRB(inFilePath,outFilePath,isLgrb)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: readDataGRB
!#endif
!
!        use pm_parallelism, only: image_type
!        use pm_kind, only: IK, LK, RK
!        implicit none
!        character(*, SK), intent(in)    :: inFilePath, outFilePath
!        integer(IK)                     :: inFileUnit, outFileUnit, igrb
!        logical(LK)     , intent(in)    :: isLgrb
!        type(image_type)                :: image
!
!        image = image_type()
!
!        if (isLgrb) then
!            GRB%count = NLGRB
!        else
!            GRB%count = NSGRB
!        end if
!
!        if (allocated(GRB%Event)) deallocate(GRB%Event); allocate(GRB%Event(GRB%count))
!        if (allocated(Trigger)) deallocate(Trigger); allocate(Trigger(GRB%count))
!
!        open( newunit = inFileUnit & ! LCOV_EXCL_LINE
!            , file = inFilePath & ! LCOV_EXCL_LINE
!            , status = "old" & ! LCOV_EXCL_LINE
!#if __INTEL_COMPILER && WINDOWS_ENABLED
!            , SHARED & ! LCOV_EXCL_LINE
!#endif
!            )
!
!        if (image%is%first) then
!        ! LCOV_EXCL_START
!            open(newunit=outFileUnit,file=outFilePath,status="replace")
!            write(outFileUnit,"(9a30)"  ) "trigger"             &
!                                        , "logPbol_1eV_20MeV"   &
!                                        , "logSbol_1eV_20MeV"   &
!                                        , "logEpk"              &
!                                        , "logEPR1024"          &
!                                        , "logEFR"              &
!                                        , "logFPR1024"          &
!                                        , "logT90"              &
!                                        , "logEffPF53"
!        end if
!        ! LCOV_EXCL_STOP
!
!        ! skip the header row in the input file
!        read(inFileUnit,*)
!
!        ! read BATSE GRB data
!        do igrb = 1, GRB%count
!
!            read(inFileUnit,*   ) Trigger(igrb)             & ! LCOV_EXCL_LINE
!                                , GRB%Event(igrb)%logPF53   &
!                                , GRB%Event(igrb)%logEpk    &
!                                , GRB%Event(igrb)%logSbol   &
!                                , GRB%Event(igrb)%logT90
!
!            ! convert all values to logarithm in base Neper
!
!            GRB%Event(igrb)%logPF53 = LOG_TEN * GRB%Event(igrb)%logPF53
!            GRB%Event(igrb)%logEpk  = LOG_TEN * GRB%Event(igrb)%logEpk
!            GRB%Event(igrb)%logSbol = LOG_TEN * GRB%Event(igrb)%logSbol
!            GRB%Event(igrb)%logT90  = LOG_TEN * GRB%Event(igrb)%logT90
!
!            ! convert photon count data to energy in units of ergs
!
!            GRB%Event(igrb)%logPbol = getLogPbol( GRB%Event(igrb)%logEpk, GRB%Event(igrb)%logPF53 )
!            if (isLgrb) then
!                GRB%Event(igrb)%logSbol = getLogPbol( GRB%Event(igrb)%logEpk, GRB%Event(igrb)%logSbol )
!            else
!                GRB%Event(igrb)%logPF53 = GRB%Event(igrb)%logPF53 - THRESH_ERFC_AMP * erfc( (GRB%Event(igrb)%logT90-THRESH_ERFC_AVG) * THRESH_ERFC_STD_INV )
!            end if
!
!            ! write the converted data to output file
!
!            if (image%is%first) then
!            ! LCOV_EXCL_START
!                write(outFileUnit,"(I30,8E30.6)") Trigger(igrb)                                     &
!                                                , GRB%Event(igrb)%logPbol                           &
!                                                , GRB%Event(igrb)%logSbol                           &
!                                                , GRB%Event(igrb)%logEpk                            &
!                                                , GRB%Event(igrb)%logEpk-GRB%Event(igrb)%logPbol    &
!                                                , GRB%Event(igrb)%logEpk-GRB%Event(igrb)%logSbol    &
!                                                , GRB%Event(igrb)%logSbol-GRB%Event(igrb)%logPbol   &
!                                                , GRB%Event(igrb)%logT90                            &
!                                                , GRB%Event(igrb)%logPF53
!
!            end if
!            ! LCOV_EXCL_STOP
!
!        end do
!
!        if (image%is%first) close(outFileUnit)
!
!        close(inFileUnit)
!
!!#if CAF_ENABLED
!!        sync images(*)
!!
!!    else
!!
!!        sync images(1)
!!        do igrb = 1, GRB%count
!!            GRB%Event(igrb)%logPbol = GRB[1]%Event(igrb)%logPbol
!!            GRB%Event(igrb)%logSbol = GRB[1]%Event(igrb)%logSbol
!!            GRB%Event(igrb)%logPF53 = GRB[1]%Event(igrb)%logPF53
!!            GRB%Event(igrb)%logEpk  = GRB[1]%Event(igrb)%logEpk
!!            GRB%Event(igrb)%logT90  = GRB[1]%Event(igrb)%logT90
!!        end do
!!
!!    end if
!!#endif
!
!    end subroutine readDataGRB
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_batse