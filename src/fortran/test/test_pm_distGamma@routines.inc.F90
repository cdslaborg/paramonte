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
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
        use pm_complexAbs, only: abs, log, operator(<), operator(<=)
#elif   !RK_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK) , parameter :: NP = 10_IK
        integer(IK) :: i
#if     RK_ENABLED || RK_ENABLED
        real(RKC)   , parameter :: TOL = epsilon(0._RKC) * 100_IK
        real(RKC)   , parameter :: invOmega = 1._RKC
        real(RKC)               :: PDF_ref(NP)
        real(RKC)               :: diff(NP)
        real(RKC)               :: PDF(NP)
        real(RKC)               :: X(NP)
        real(RKC)               :: Kappa(NP)
        real(RKC)               :: InvSigma(NP)
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK
        call setUnifRand(X, epsilon(0._RKC), 1000._RKC)
        call setUnifRand(Kappa, epsilon(0._RKC), 10._RKC)
        call setUnifRand(InvSigma, epsilon(0._RKC), 10._RKC)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PDF_ref = getGenGammaLogPDF(X)
#if     getGammaLogPDF_ENABLED
        PDF = getGammaLogPDF(X)
#elif   setGammaLogPDF_ENABLED
        call setGammaLogPDF(PDF, X)
#else
#error  "Unrecognized interface."
#endif
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PDF_ref = getGenGammaLogPDF(X, Kappa)
#if     getGammaLogPDF_ENABLED
        PDF = getGammaLogPDF(X, Kappa)
#elif   setGammaLogPDF_ENABLED
        call setGammaLogPDF(PDF, X, getGammaLogPDFNF(Kappa), Kappa)
#else
#error  "Unrecognized interface."
#endif
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PDF_ref = getGenGammaLogPDF(X, Kappa, InvSigma = InvSigma)
#if     getGammaLogPDF_ENABLED
        PDF = getGammaLogPDF(X, Kappa, InvSigma)
#elif   setGammaLogPDF_ENABLED
        call setGammaLogPDF(PDF, X, getGammaLogPDFNF(Kappa, InvSigma), Kappa, InvSigma)
#else
#error  "Unrecognized interface."
#endif
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer(IK), intent(in) :: line
            diff = abs(PDF - PDF_ref)
            do i = 1, NP
                call test%assert(assertion, desc = "The PDF of the Gamma distribution must be computed correctly.", line = line)
                assertion = assertion .and. diff(i) <= TOL
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "InvSigma   ", InvSigma(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "Kappa      ", Kappa(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "PDF_ref    ", PDF_ref(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "PDF        ", PDF(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "X          ", X(i)
                    write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff(i)
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
