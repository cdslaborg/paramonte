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
!>  This include file contains procedure implementations of the tests of [pm_distGamma](@ref pm_distGamma).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
        use pm_complexAbs, only: abs, log, operator(<), operator(<=)
#elif   !RK_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK) :: i
        integer(IK) , parameter :: NP = 10_IK
        real(TKC)   , parameter :: TOL = epsilon(0._TKC) * 100
        real(TKC)   , parameter :: invOmega = 1._TKC
        real(TKC)               :: pdf_ref(NP)
        real(TKC)               :: diff(NP)
        real(TKC)               :: pdf(NP)
        real(TKC)               :: point(NP)
        real(TKC)               :: kappa(NP)
        real(TKC)               :: invSigma(NP)

        assertion = .true._LK
        call setUnifRand(point, epsilon(0._TKC), 1000._TKC)
        call setUnifRand(kappa, epsilon(0._TKC), 10._TKC)
        call setUnifRand(invSigma, epsilon(0._TKC), 10._TKC)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pdf_ref = getGenGammaLogPDF(point)
#if     getGammaLogPDF_ENABLED
        pdf = getGammaLogPDF(point)
#elif   setGammaLogPDF_ENABLED
        call setGammaLogPDF(pdf, point)
#else
#error  "Unrecognized interface."
#endif
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pdf_ref = getGenGammaLogPDF(point, kappa)
#if     getGammaLogPDF_ENABLED
        pdf = getGammaLogPDF(point, kappa)
#elif   setGammaLogPDF_ENABLED
        call setGammaLogPDF(pdf, point, getGammaLogPDFNF(kappa), kappa)
#else
#error  "Unrecognized interface."
#endif
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pdf_ref = getGenGammaLogPDF(point, kappa, invSigma = invSigma)
#if     getGammaLogPDF_ENABLED
        pdf = getGammaLogPDF(point, kappa, invSigma)
#elif   setGammaLogPDF_ENABLED
        call setGammaLogPDF(pdf, point, getGammaLogPDFNF(kappa, invSigma), kappa, invSigma)
#else
#error  "Unrecognized interface."
#endif
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer(IK), intent(in) :: line
            diff = abs(pdf - pdf_ref)
            do i = 1, NP
                call test%assert(assertion, desc = "The pdf of the Gamma distribution must be computed correctly.", line = line)
                assertion = isClose(pdf_ref(i), pdf(i), reltol = TOL) .and. assertion
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "invSigma   ", invSigma(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "kappa      ", kappa(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "pdf_ref    ", pdf_ref(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "pdf        ", pdf(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "point      ", point(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                    write(test%disp%unit,"(*(g0,:,', '))") "diff(i)<TOL", diff(i)<TOL
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
