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
!>  This include file contains procedure implementations of the tests of [pm_cosmology](@ref pm_cosmology).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !%%%%%%%%%%%%%%%%%%%%%%%%
#if         getSizeUnivNormed_ENABLED
            !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKC)   , parameter :: TOL = epsilon(0._RKC)*100
        real(RKC)               :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKC)               :: hubbleParamNormedSq_ref
        real(RKC)               :: hubbleParamNormedSq
        real(RKC)               :: diff

        assertion = .true._LK
        do i = 1_IK, 500_IK

            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            zplus1 = getUnifRand(1._RKC, 100._RKC)
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1)
            call report(int(__LINE__, IK))

            omegaM = getUnifRand(0._RKC, 1._RKC)
            omegaL = 1._RKC - omegaM
            zplus1 = getUnifRand(1._RKC, 100._RKC)
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL)
            call report(int(__LINE__, IK))

            omegaM = getUnifRand(0._RKC, 1._RKC)
            omegaL = merge(getUnifRand(0._RKC, 1._RKC - omegaM), 0._RKC, 1._RKC - omegaM > 0._RKC)
            omegaR = 1._RKC - omegaM - omegaL
            zplus1 = getUnifRand(1._RKC, 100._RKC)
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)
            call report(int(__LINE__, IK))

            omegaM = getUnifRand(0._RKC, 1._RKC)
            omegaL = merge(getUnifRand(0._RKC, 1._RKC - omegaM), 0._RKC, 1._RKC - omegaM > 0._RKC)
            omegaR = merge(getUnifRand(0._RKC, 1._RKC - omegaM - omegaL), 0._RKC, 1._RKC - omegaM - omegaL > 0._RKC)
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)
            call report(int(__LINE__, IK))

        end do

    contains

        function getHubbleParamNormedSq_ref(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq_ref)
            real(RKC), intent(in) :: zplus1, omegaM, omegaL, omegaR, omegaK
            real(RKC) :: hubbleParamNormedSq_ref
            hubbleParamNormedSq_ref = omegaL + zplus1**2 * (omegaK + zplus1 * (omegaM + zplus1 * omegaR))
        end function

        subroutine report(line)
            integer(IK), intent(in) :: line
            hubbleParamNormedSq_ref = getHubbleParamNormedSq_ref(zplus1, omegaM, omegaL, omegaR, omegaK)
            diff = abs(hubbleParamNormedSq - hubbleParamNormedSq_ref) / hubbleParamNormedSq_ref
            assertion = assertion .and. logical(diff < TOL, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "hubbleParamNormedSq_ref", hubbleParamNormedSq_ref
                write(test%disp%unit,"(*(g0,:,', '))") "hubbleParamNormedSq    ", hubbleParamNormedSq
                write(test%disp%unit,"(*(g0,:,', '))") "zplus1                 ", zplus1
                write(test%disp%unit,"(*(g0,:,', '))") "omegaM                 ", omegaM
                write(test%disp%unit,"(*(g0,:,', '))") "omegaL                 ", omegaL
                write(test%disp%unit,"(*(g0,:,', '))") "omegaR                 ", omegaR
                write(test%disp%unit,"(*(g0,:,', '))") "omegaK                 ", omegaK
                write(test%disp%unit,"(*(g0,:,', '))") "diff                   ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL                    ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "The dimensionless Hubble Parameter must be computed correctly for the specified cosmological parameters.", line = line)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisLookbackNormed_ENABLED || getDisComNormed_ENABLED || getDisComTransNormed_ENABLED || \
        getDisAngNormed_ENABLED || getDisLumNormed_ENABLED || getDisComTransNormedWU10_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK)             :: neval, err
        real(RKC)   , parameter :: reltol = epsilon(0._RKC) * 10000!sqrt(epsilon(0._RKC))
        real(RKC)   , parameter :: MPC2LY_RKC = real(MPC2LY, RKC)
        real(RKC)   , parameter :: LIGHT_SPEED_RKC = real(LIGHT_SPEED, RKC)
        real(RKC)   , parameter :: HUBBLE_CONST_RKC = real(HUBBLE_CONST, RKC)
        real(RKC)   , parameter :: HUBBLE_DISTANCE_MPC_RKC = real(HUBBLE_DISTANCE_MPC, RKC)
        real(RKC)               :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKC)               :: distance_ref
        real(RKC)               :: distance
        real(RKC)               :: diff
        real(RKC)               :: tol

        tol = epsilon(0._RKC) * 10000
        assertion = .true._LK

        call runTestsWith()
        call runTestsWith(neval)
        call runTestsWith(err = err)
        call runTestsWith(neval, err)

    contains

        subroutine runTestsWith(neval, err)
            integer(IK), intent(out), optional :: neval, err

#if         getSizeUnivNormed_ENABLED

            zplus1 = 1._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            distance_ref = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) ! 0.955239468086737240441960122539593084_RKC

            distance = getSizeUnivNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = huge(0._RKC)
            distance_ref = 0._RKC

            distance = getSizeUnivNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 3.5_RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = 0.560176566313766147208125382717677278E-1_RKC

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = +4._RKC
            omegaM = +0.2_RKC
            omegaL = -0.5_RKC
            omegaR = +0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR
            distance_ref = 0.418795901922915638929080324595150918E-1_RKC

            distance = getSizeUnivNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Universe Size must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

#elif       getDisLookbackNormed_ENABLED

            zplus1 = 1001._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            distance_ref = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) ! 0.955239468086737240441960122539593084_RKC

            distance = getDisLookbackNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 0._RKC
            distance_ref = 0._RKC

            distance = getDisLookbackNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 6._RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = 0.586437343702725508511169269833928737_RKC

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = +4._RKC
            omegaM = +0.2_RKC
            omegaL = -0.5_RKC
            omegaR = +0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR
            distance_ref = 0.501459151132179770407237268671874632_RKC

            distance = getDisLookbackNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Lookback Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

#elif       getDisComNormed_ENABLED

            zplus1 = 2._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            distance_ref = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) ! 3396.03639669472458732702995144103428_RKC

            distance = getDisComNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 0._RKC
            distance_ref = 0._RKC

            distance = getDisComNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 6._RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = 4608.74661107738061825469278425187163_RKC / HUBBLE_DISTANCE_MPC_RKC

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 6._RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = 4608.74661107738061825469278425187163_RKC / HUBBLE_DISTANCE_MPC_RKC

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = +4._RKC
            omegaM = +0.2_RKC
            omegaL = -0.5_RKC
            omegaR = +0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR
            distance_ref = 3656.22113700199274111901199467371355_RKC / HUBBLE_DISTANCE_MPC_RKC

            distance = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

#elif       getDisComTransNormed_ENABLED

            zplus1 = 2._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)

            distance = getDisComTransNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Transverse Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Transverse Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Transverse Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Transverse Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = +4._RKC
            omegaM = +0.2_RKC
            omegaL = -0.5_RKC
            omegaR = +0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)

            distance = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Transverse Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 6._RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)

            distance = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Transverse Comoving Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

#elif       getDisAngNormed_ENABLED

            zplus1 = 2._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) / zplus1

            distance = getDisAngNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Angular Diameter Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisAngNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Angular Diameter Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisAngNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Angular Diameter Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisAngNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Angular Diameter Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = +4._RKC
            omegaM = +0.2_RKC
            omegaL = -0.5_RKC
            omegaR = +0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) / zplus1

            distance = getDisAngNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Angular Diameter Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 6._RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) / zplus1

            distance = getDisAngNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Angular Diameter Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

#elif       getDisLumNormed_ENABLED

            zplus1 = 2._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) * zplus1

            distance = getDisLumNormed(zplus1, reltol, neval, err)
            call report()
            call test%assert(assertion, SK_"The Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLumNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLumNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            distance = getDisLumNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = +4._RKC
            omegaM = +0.2_RKC
            omegaL = -0.5_RKC
            omegaR = +0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) * zplus1

            distance = getDisLumNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 6._RKC
            omegaM = 0.2_RKC
            omegaL = 0.5_RKC
            omegaR = 0.5_RKC
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            distance_ref = getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) * zplus1

            distance = getDisLumNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#elif       getDisComTransNormedWU10_ENABLED

            tol = 0.002_RKC
            block
                integer :: i
                do i = 1, 5
                    zplus1 = 1.1_RKC * i
                    omegaM = real(OMEGA_M, RKC)
                    omegaL = real(OMEGA_L, RKC)
                    omegaR = 0._RKC
                    omegaK = 0._RKC
                    distance_ref = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)

                    distance = getDisComTransNormedWU10(zplus1)
                    call report()
                    call test%assert(assertion, SK_"The WU10 Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                    distance = getDisComTransNormedWU10(zplus1, omegaM, omegaL)
                    call report(omegaM, omegaL)
                    call test%assert(assertion, SK_"The WU10 Luminosity Distance must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))
                end do
            end block

#else
#error  "Unrecognized interface."
#endif
        end subroutine

#if     getDisComTransNormed_ENABLED || getDisAngNormed_ENABLED || getDisLumNormed_ENABLED
        function getDisComTransNormed_def(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err) result(disComTransNormed)
            real(RKC), intent(in), optional :: omegaM, omegaL, omegaR, omegaK, reltol
            integer(IK), intent(out) :: neval, err
            real(RKC), intent(in) :: zplus1
            real(RKC):: disComTransNormed
            real(RKC):: sqrtOmegaK
            if (present(omegaM) .and. present(omegaL) .and. present(omegaR) .and. present(omegaK)) then
                disComTransNormed = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
            elseif (present(omegaM) .and. present(omegaL) .and. present(omegaR)) then
                disComTransNormed = getDisComNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
            elseif (present(omegaM) .and. present(omegaL)) then
                disComTransNormed = getDisComNormed(zplus1, omegaM, omegaL, reltol, neval, err)
            else
                disComTransNormed = getDisComNormed(zplus1, reltol, neval, err)
            end if
            if (present(omegaK)) then
                if (omegaK > epsilon(0._RKC)) then
                    sqrtOmegaK = sqrt(omegaK)
                    disComTransNormed = sinh(disComTransNormed * sqrtOmegaK) / sqrtOmegaK
                elseif (omegaK < -epsilon(0._RKC)) then
                    sqrtOmegaK = sqrt(-omegaK)
                    disComTransNormed = sin(disComTransNormed * sqrtOmegaK) / sqrtOmegaK
                end if
            end if
        end function
#endif
        subroutine report(omegaM, omegaL, omegaR, omegaK)
            real(RKC), intent(in), optional :: omegaM, omegaL, omegaR, omegaK
            if (abs(distance_ref) <= epsilon(0._RKC)) then
                diff = abs(distance - distance_ref)
            else
                diff = abs(distance - distance_ref) / distance_ref
            end if
            assertion = assertion .and. logical(diff < TOL, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "distance_ref", distance_ref
                write(test%disp%unit,"(*(g0,:,', '))") "distance", distance
                write(test%disp%unit,"(*(g0,:,', '))") "zplus1", zplus1
                if (present(omegaM)) write(test%disp%unit,"(*(g0,:,', '))") "omegaM", omegaM
                if (present(omegaL)) write(test%disp%unit,"(*(g0,:,', '))") "omegaL", omegaL
                if (present(omegaR)) write(test%disp%unit,"(*(g0,:,', '))") "omegaR", omegaR
                if (present(omegaK)) write(test%disp%unit,"(*(g0,:,', '))") "omegaK", omegaK
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getHubbleParamNormedSq_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKC)   , parameter :: TOL = epsilon(0._RKC)*100
        real(RKC)               :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKC)               :: hubbleParamNormedSq_ref
        real(RKC)               :: hubbleParamNormedSq
        real(RKC)               :: diff
        assertion = .true._LK

        do i = 1_IK, 500_IK

            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            zplus1 = getUnifRand(1._RKC, 100._RKC)
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1)
            call report(int(__LINE__, IK))

            omegaM = getUnifRand(0._RKC, 1._RKC)
            omegaL = 1._RKC - omegaM
            zplus1 = getUnifRand(1._RKC, 100._RKC)
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL)
            call report(int(__LINE__, IK))

            omegaM = getUnifRand(0._RKC, 1._RKC)
            omegaL = merge(getUnifRand(0._RKC, 1._RKC - omegaM), 0._RKC, 1._RKC - omegaM > 0._RKC)
            omegaR = 1._RKC - omegaM - omegaL
            zplus1 = getUnifRand(1._RKC, 100._RKC)
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)
            call report(int(__LINE__, IK))

            omegaM = getUnifRand(0._RKC, 1._RKC)
            omegaL = merge(getUnifRand(0._RKC, 1._RKC - omegaM), 0._RKC, 1._RKC - omegaM > 0._RKC)
            omegaR = merge(getUnifRand(0._RKC, 1._RKC - omegaM - omegaL), 0._RKC, 1._RKC - omegaM - omegaL > 0._RKC)
            omegaK = 1._RKC - omegaM - omegaL - omegaR
            hubbleParamNormedSq = getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)
            call report(int(__LINE__, IK))

        end do

    contains

        function getHubbleParamNormedSq_ref(zplus1, omegaM, omegaL, omegaR, omegaK) result(hubbleParamNormedSq_ref)
            real(RKC), intent(in) :: zplus1, omegaM, omegaL, omegaR, omegaK
            real(RKC) :: hubbleParamNormedSq_ref
            hubbleParamNormedSq_ref = omegaL + zplus1**2 * (omegaK + zplus1 * (omegaM + zplus1 * omegaR))
        end function

        subroutine report(line)
            integer(IK), intent(in) :: line
            hubbleParamNormedSq_ref = getHubbleParamNormedSq_ref(zplus1, omegaM, omegaL, omegaR, omegaK)
            diff = abs(hubbleParamNormedSq - hubbleParamNormedSq_ref) / hubbleParamNormedSq_ref
            assertion = assertion .and. logical(diff < TOL, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "hubbleParamNormedSq_ref", hubbleParamNormedSq_ref
                write(test%disp%unit,"(*(g0,:,', '))") "hubbleParamNormedSq    ", hubbleParamNormedSq    
                write(test%disp%unit,"(*(g0,:,', '))") "zplus1                 ", zplus1
                write(test%disp%unit,"(*(g0,:,', '))") "omegaM                 ", omegaM
                write(test%disp%unit,"(*(g0,:,', '))") "omegaL                 ", omegaL
                write(test%disp%unit,"(*(g0,:,', '))") "omegaR                 ", omegaR
                write(test%disp%unit,"(*(g0,:,', '))") "omegaK                 ", omegaK
                write(test%disp%unit,"(*(g0,:,', '))") "diff                   ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL                    ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "The dimensionless Hubble Parameter must be computed correctly for the specified cosmological parameters.", line = line)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getVolComDiffNormed_ENABLED || setVolComDiffNormed_ENABLED || getVolComNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK)             :: neval, err
        real(RKC)   , parameter :: reltol = epsilon(0._RKC) * 10000!sqrt(epsilon(0._RKC))
        real(RKC)   , parameter :: MPC2LY_RKC = real(MPC2LY, RKC)
        real(RKC)   , parameter :: LIGHT_SPEED_RKC = real(LIGHT_SPEED, RKC)
        real(RKC)   , parameter :: HUBBLE_CONST_RKC = real(HUBBLE_CONST, RKC)
        real(RKC)   , parameter :: HUBBLE_DISTANCE_MPC_RKC = real(HUBBLE_DISTANCE_MPC, RKC)
        real(RKC)               :: zplus1, omegaM, omegaL, omegaR, omegaK
        real(RKC)               :: volComNormed
        real(RKC)               :: volComNormed_ref
        real(RKC)               :: diff
        real(RKC)               :: tol

        tol = epsilon(0._RKC) * 10000
        assertion = .true._LK

        call runTestsWith()
        call runTestsWith(neval)
        call runTestsWith(err = err)
        call runTestsWith(neval, err)

    contains

        subroutine runTestsWith(neval, err)
            use pm_distUnif, only: setUnifRand
            integer(IK), intent(out), optional :: neval, err

#if         getVolComDiffNormed_ENABLED

            integer(IK) :: i
            do i = 1, 10

                omegaM = real(OMEGA_M, RKC)
                omegaL = real(OMEGA_L, RKC)
                omegaR = real(OMEGA_R, RKC)
                omegaK = real(OMEGA_K, RKC)
                call setUnifRand(zplus1, 1._RKC, 10._RKC)
                volComNormed = getVolComDiffNormed(zplus1, reltol, neval, err)
                call setVolComDiffNormed(volComNormed_ref, getDisComTransNormed(zplus1, reltol)**2, sqrt(getHubbleParamNormedSq(zplus1)))
                call report(omegaM, omegaL, omegaR, omegaK)
                call test%assert(assertion, SK_"@getVolComDiffNormed(): The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                call setUnifRand(omegaM, 0._RKC, 1._RKC)
                omegaL = 1._RKC - omegaM
                omegaR = real(OMEGA_R, RKC)
                omegaK = real(OMEGA_K, RKC)
                call setUnifRand(zplus1, 1._RKC, 10._RKC)
                volComNormed = getVolComDiffNormed(zplus1, omegaM, omegaL, reltol, neval, err)
                call setVolComDiffNormed(volComNormed_ref, getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL)))
                call report(omegaM, omegaL, omegaR, omegaK)
                call test%assert(assertion, SK_"@getVolComDiffNormed(): The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                call setUnifRand(omegaM, 0._RKC, 1._RKC)
                call setUnifRand(omegaL, 0._RKC, 1._RKC - omegaM)
                omegaR = 1._RKC - omegaM - omegaL
                omegaK = real(OMEGA_K, RKC)
                call setUnifRand(zplus1, 1._RKC, 10._RKC)
                volComNormed = getVolComDiffNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
                call setVolComDiffNormed(volComNormed_ref, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)))
                call report(omegaM, omegaL, omegaR, omegaK)
                call test%assert(assertion, SK_"@getVolComDiffNormed(): The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                call setUnifRand(omegaM, 0._RKC, 1._RKC)
                call setUnifRand(omegaL, 0._RKC, omegaM)
                call setUnifRand(omegaR, 0._RKC, omegaM + omegaL)
                omegaK = 1._RKC - omegaM - omegaL - omegaR
                call setUnifRand(zplus1, 1._RKC, 10._RKC)
                volComNormed = getVolComDiffNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)
                call setVolComDiffNormed(volComNormed_ref, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
                call report(omegaM, omegaL, omegaR, omegaK)
                call test%assert(assertion, SK_"@getVolComDiffNormed(): The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            end do

#elif       setVolComDiffNormed_ENABLED

            omegaM = 0.2_RKC
            omegaL = 0.4_RKC
            omegaR = 0.0_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR

            zplus1 = 1._RKC
            volComNormed_ref = 0._RKC
            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            zplus1 = 2.5_RKC
            volComNormed_ref = 0.424491473544200570744727780980436799_RKC
            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            zplus1 = 6._RKC
            volComNormed_ref = 0.601981403654138797418469447630789515_RKC
            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 1._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            volComNormed_ref = 0._RKC

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1)))
            call report()
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL)))
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)))
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 2._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            call setVolComDiffNormed_def(volComNormed_ref, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1)))
            call report()
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL)))
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)))
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            zplus1 = 3._RKC
            omegaM = real(OMEGA_M, RKC)
            omegaL = real(OMEGA_L, RKC)
            omegaR = real(OMEGA_R, RKC)
            omegaK = real(OMEGA_K, RKC)
            call setVolComDiffNormed_def(volComNormed_ref, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1)))
            call report()
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL)))
            call report(omegaM, omegaL)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)))
            call report(omegaM, omegaL, omegaR)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            call setVolComDiffNormed(volComNormed, getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err)**2, sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

#elif       getVolComNormed_ENABLED

            omegaM = 0.2_RKC
            omegaL = 0.4_RKC
            omegaR = 0.0_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR

            zplus1 = 1._RKC
            volComNormed_ref = 0._RKC
            volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            zplus1 = 2.5_RKC
            volComNormed_ref = 3.99648746549472192582313390774500827_RKC
            volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            zplus1 = 6._RKC
            volComNormed_ref = 29.0231597572974447624597993854436293_RKC
            volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            omegaM = 0.2_RKC
            omegaL = 0.4_RKC
            omegaR = 0.5_RKC
            omegaK = +1._RKC - omegaM - omegaL - omegaR

            zplus1 = 1._RKC
            volComNormed_ref = 0._RKC
            volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            zplus1 = 2.5_RKC
            volComNormed_ref = 1.50955035302854263049196839502469953_RKC
            volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            zplus1 = 6._RKC
            volComNormed_ref = 4.45878362364779446349477242289460663_RKC
            volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
            call report(omegaM, omegaL, omegaR, omegaK)
            call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            block
                integer(IK) :: i
                do i = 1, 3
                    zplus1 = real(i, RKC)
                    omegaM = real(OMEGA_M, RKC)
                    omegaL = real(OMEGA_L, RKC)
                    omegaR = real(OMEGA_R, RKC)
                    omegaK = real(OMEGA_K, RKC)
                    volComNormed_ref = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))

                    volComNormed = getVolComNormed(getDisComTransNormed(zplus1, reltol, neval, err))
                    call report()
                    call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                    volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err))
                    call report(omegaM, omegaL)
                    call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                    volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err))
                    call report(omegaM, omegaL, omegaR)
                    call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))

                    volComNormed = getVolComNormed(getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrt(abs(omegaK)), reltol, neval, err), omegaK, sqrt(abs(omegaK)))
                    call report(omegaM, omegaL, omegaR, omegaK)
                    call test%assert(assertion, SK_"The Comoving Volume Element must be computed correctly for the specified cosmological parameters with present(neval), present(err) = "//getStr([present(neval), present(err)]), int(__LINE__, IK))
                end do
            end block
#else
#error  "Unrecognized interface."
#endif
        end subroutine

#if     setVolComDiffNormed_ENABLED
        pure elemental subroutine setVolComDiffNormed_def(volComNormed, disComTransNormedSq, hubbleParamNormed)
            real(RKC), intent(in)   :: disComTransNormedSq, hubbleParamNormed
            real(RKC), intent(out)  :: volComNormed
            volComNormed = disComTransNormedSq / hubbleParamNormed
        end subroutine
#endif
        subroutine report(omegaM, omegaL, omegaR, omegaK)
            real(RKC), intent(in), optional :: omegaM, omegaL, omegaR, omegaK
            if (abs(volComNormed_ref) <= epsilon(0._RKC)) then
                diff = abs(volComNormed - volComNormed_ref)
            else
                diff = abs(volComNormed - volComNormed_ref) / volComNormed_ref
            end if
            assertion = assertion .and. logical(diff < TOL, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "volComNormed_ref", volComNormed_ref
                write(test%disp%unit,"(*(g0,:,', '))") "volComNormed", volComNormed
                write(test%disp%unit,"(*(g0,:,', '))") "zplus1", zplus1
                if (present(omegaM)) write(test%disp%unit,"(*(g0,:,', '))") "omegaM", omegaM
                if (present(omegaL)) write(test%disp%unit,"(*(g0,:,', '))") "omegaL", omegaL
                if (present(omegaR)) write(test%disp%unit,"(*(g0,:,', '))") "omegaR", omegaR
                if (present(omegaK)) write(test%disp%unit,"(*(g0,:,', '))") "omegaK", omegaK
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif