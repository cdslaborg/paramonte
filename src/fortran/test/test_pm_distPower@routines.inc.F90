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
!>  This include file contains procedure implementations of the tests of [pm_distPower](@ref pm_distPower).
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
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getPowerLogPDF_ENABLED || setPowerLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = sqrt(epsilon(0._RKG))
        real(RKG)                   :: alpha, logPDFNF, logMinX, logMaxX
        real(RKG)                   :: output, output_ref, diff

        assertion = .true._LK

        do i = 1_IK, 50_IK
            logMinX = getUnifRand(-5._RKG, +4._RKG)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMinX)
            call runTestsWith()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(logMinX)

            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(log(huge(0._RKG)))
            real(RKG), intent(in), optional :: logMinX

            alpha = getUnifRand(0.5_RKG, +3._RKG)
            if (present(logMinX)) then
                logPDFNF = getPowerLogPDFNF(alpha, logMinX, logMaxX)
            else
                logPDFNF = getPowerLogPDFNF(alpha, logMaxX)
            end if

#if         getPowerLogPDF_ENABLED
            block
                real(RKG) :: logx
                logx = getUnifRand(getOption(-SQRT_LOG_HUGE, logMinX), logMaxX)
                call setPowerLogPDF(output_ref, logx, alpha, logPDFNF)
                if (present(logMinX)) then
                    output = getPowerLogPDF(logx, alpha, logMinX, logMaxX)
                else
                    output = getPowerLogPDF(logx, alpha, logMaxX)
                end if
                call compare(logMinX, logx)
            end block
#elif       setPowerLogPDF_ENABLED
            output_ref = 1._RKG
            block
                use pm_quadPack, only: isFailedQuad, getQuadErr, GK31, weps
                integer(IK) :: err, neval, nint, sindex(2000)
                real(RKG) :: lb, abserr, sinfo(4,2000)
                character(255, SK) :: msg
                msg = SK_" "
                lb = 0._RKG
                if (present(logMinX)) lb = exp(logMinX)
                if (alpha > 1._RKG) then
                    assertion = assertion .and. .not. isFailedQuad(getPowerPDF, lb, exp(logMaxX), output, reltol = TOL, msg = msg)
                else
                    err = getQuadErr(getPowerPDF, lb, exp(logMaxX), 0._RKG, TOL, GK31, weps, output, abserr, sinfo, sindex, neval, nint) ! extends QAGS/QAGI routines of QuadPack.
                    assertion = assertion .and. err == 0_IK
                end if
                call report(logMinX)
                call test%assert(assertion, SK_"@setPowerLogPDF(): The integral of the Power distribution must be computed without failure. msg: "//trim(msg), int(__LINE__, IK))
            end block
            call compare()
            call report(logMinX)
            call test%assert(assertion, SK_"@setPowerLogPDF(): The PDF of the Power distribution must be computed correctly.", int(__LINE__, IK))
#else
#error      "Unrecognized interface."
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getPowerLogPDF_ENABLED
        subroutine compare(logMinX, logx)
            real(RKG), intent(in), optional :: logMinX, logx
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
            call report(logMinX, logx)
            call test%assert(assertion, SK_"@getPowerLogPDF(): The PDF of the Power distribution must be computed correctly.", int(__LINE__, IK))
        end subroutine
#elif   setPowerLogPDF_ENABLED
        subroutine compare()
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
        end subroutine
        function getPowerPDF(x) result(pdf)
            real(RKG), intent(in) :: x
            real(RKG) :: pdf
            call setPowerLogPDF(pdf, log(x), alpha, logPDFNF)
            pdf = exp(pdf)
        end function
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(logMinX, logx)
            real(RKG), intent(in), optional :: logMinX, logx
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(logx)", present(logx)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMinX)", present(logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, logx)", getOption(HUGE_RKG, logx)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMinX)", getOption(-HUGE_RKG, logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "logPDFNF", logPDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMaxX", logMaxX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref", output_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogCDF_ENABLED || setPowerLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = sqrt(epsilon(0._RKG))
        real(RKG)                   :: alpha, logCDFNF, logMinX, logMaxX
        real(RKG)                   :: output, output_ref, diff
        logical(LK)                 :: logMinXPresent

        assertion = .true._LK

        do i = 1_IK, 300_IK
            logMinX = getUnifRand(-5._RKG, +4._RKG)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMinX = logMinX)
            call runTestsWith()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(logMinX)

            use pm_val2str, only: getStr
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(LOG_HUGE)
            real(RKG), intent(in), optional :: logMinX
            real(RKG) :: logx

            logMinXPresent = present(logMinX)
            alpha = getUnifRand(0.5_RKG, +3._RKG)
            if (logMinXPresent) then
                logCDFNF = getPowerLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = getPowerLogCDFNF(alpha, logMaxX)
            end if

#if         getPowerLogCDF_ENABLED

            logx = getUnifRand(getOption(-SQRT_LOG_HUGE, logMinX), logMaxX)

            if (logMinXPresent) then
                call setPowerLogCDF(output_ref, logx, alpha, logMinX, logCDFNF)
                output = getPowerLogCDF(logx, alpha, logMinX, logMaxX)
            else
                call setPowerLogCDF(output_ref, logx, alpha, logCDFNF)
                output = getPowerLogCDF(logx, alpha, logMaxX)
            end if
            call compare(__LINE__, logMinX, logx)

#elif       setPowerLogCDF_ENABLED

            output_ref = -LOG_HUGE
            logx = getOption(-LOG_HUGE, logMinX)

            if (logMinXPresent) then
                call setPowerLogCDF(output, logx, alpha, logMinX, logCDFNF)
            else
                call setPowerLogCDF(output, logx, alpha, logCDFNF)
            end if
            call report(logMinX, logx)
            call test%assert(assertion, SK_"@setPowerLogCDF(): The CDF of the Power distribution at the lower limit of the support must be zero.", int(__LINE__, IK))

            output_ref = 0._RKG
            logx = logMaxX

            if (logMinXPresent) then
                call setPowerLogCDF(output, logx, alpha, logMinX, logCDFNF)
            else
                call setPowerLogCDF(output, logx, alpha, logCDFNF)
            end if
            call report(logMinX, logx)
            call test%assert(assertion, SK_"@setPowerLogCDF(): The CDF of the Power distribution at the upper limit of the support must be unity.", int(__LINE__, IK))

            ! Compare the CDF difference at random points with the corresponding integral of the PDF.

            block
                use pm_quadPack, only: isFailedQuad, getQuadErr, GK31, weps
                integer(IK) :: err, neval, nint, sindex(2000)
                real(RKG) :: lb, abserr, sinfo(4,2000)
                character(255, SK) :: msg
                msg = SK_" "
                lb = 0._RKG
                if (logMinXPresent) lb = exp(logMinX)
                logx = log(getUnifRand((exp(logMaxX) + lb) * .5_RKG, exp(logMaxX)))
                if (alpha > 1._RKG) then
                    assertion = assertion .and. .not. isFailedQuad(getPowerPDF, lb, exp(logx), output_ref, reltol = TOL, msg = msg)
                else
                    err = getQuadErr(getPowerPDF, lb, exp(logx), 0._RKG, TOL, GK31, weps, output_ref, abserr, sinfo, sindex, neval, nint) ! extends QAGS/QAGI routines of QuadPack.
                    assertion = assertion .and. err == 0_IK
                end if
                call report(logMinX, logx)
                call test%assert(assertion, SK_"@setPowerLogCDF(): The integral of the Power distribution must be computed without failure. msg: "//trim(msg), int(__LINE__, IK))
                if (logMinXPresent) then
                    call setPowerLogCDF(output, logx, alpha, logMinX, logCDFNF)
                else
                    call setPowerLogCDF(output, logx, alpha, logCDFNF)
                end if
                output = exp(output)
                call compare()
                call report(logMinX, logx)
                call test%assert(assertion, SK_"@setPowerLogCDF(): The CDF of the Power distribution must be computed correctly.", int(__LINE__, IK))
            end block

#else
#error      "Unrecognized interface."
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getPowerLogCDF_ENABLED
        subroutine compare(line, logMinX, logx)
            integer :: line
            real(RKG), intent(in), optional :: logMinX, logx
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
            call report(logMinX, logx)
            call test%assert(assertion, SK_"@getPowerLogCDF(): The CDF of the Power distribution must be computed correctly.", int(line, IK))
        end subroutine
#elif   setPowerLogCDF_ENABLED
        subroutine compare()
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
        end subroutine
        function getPowerPDF(x) result(pdf)
            use pm_distPower, only: getPowerLogPDF
            real(RKG), intent(in) :: x
            real(RKG) :: pdf
            if (logMinXPresent) then
                pdf = exp(getPowerLogPDF(log(x), alpha, logMinX, logMaxX))
            else
                pdf = exp(getPowerLogPDF(log(x), alpha, logMaxX))
            end if
        end function
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(logMinX, logx)
            real(RKG), intent(in), optional :: logMinX, logx
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(logx)", present(logx)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMinX)", present(logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, logx)", getOption(HUGE_RKG, logx)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMinX)", getOption(-HUGE_RKG, logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "logCDFNF", logCDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMaxX", logMaxX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref", output_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogQuan_ENABLED || setPowerLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = sqrt(epsilon(0._RKG))! * 1000
        real(RKG)                   :: alpha, logCDFNF, logMinX, logMaxX
        real(RKG)                   :: logx, output, output_ref, diff
        logical(LK)                 :: logMinXPresent

        assertion = .true._LK

        do i = 1_IK, 500_IK
            logMinX = getUnifRand(-LOG_HUGE, LOG_HUGE)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMinX = logMinX)
            call runTestsWith()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(logMinX)

            use pm_val2str, only: getStr
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(LOG_HUGE)
            real(RKG), intent(in), optional :: logMinX

            logMinXPresent = present(logMinX)
            alpha = getUnifRand(0.5_RKG, +3._RKG)
            if (logMinXPresent) then
                logCDFNF = getPowerLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = getPowerLogCDFNF(alpha, logMaxX)
            end if

            if (i > 1_IK) then
                logx = getUnifRand(getOption(logMaxX - 4._RKG, logMinX), logMaxX)
            else ! test for logCDF = 0.
                logx = logMaxX
                output_ref = 0._RKG ! logCDF
            end if
            if (logMinXPresent) then
                if (i > 1_IK) call setPowerLogCDF(output_ref, logx, alpha, logMinX, logCDFNF)
#if             getPowerLogQuan_ENABLED
                output = getPowerLogQuan(output_ref, alpha, logMinX, logMaxX)
#elif           setPowerLogQuan_ENABLED
                call setPowerLogQuan(output, output_ref, alpha, logMinX, logCDFNF)
#else
#error          "Unrecognized interface."
#endif
            else
                if (i > 1_IK) call setPowerLogCDF(output_ref, logx, alpha, logCDFNF)
#if             getPowerLogQuan_ENABLED
                output = getPowerLogQuan(output_ref, alpha, logMaxX)
#elif           setPowerLogQuan_ENABLED
                call setPowerLogQuan(output, output_ref, alpha, logCDFNF)
#endif
            end if
            output_ref = logx
            call report(__LINE__, logMinX, logx)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, logMinX, logx)
            integer, intent(in) :: line
            real(RKG), intent(in), optional :: logMinX, logx
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(logx)", present(logx)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMinX)", present(logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, logx)", getOption(HUGE_RKG, logx)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMinX)", getOption(-HUGE_RKG, logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "logCDFNF", logCDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMaxX", logMaxX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref", output_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@getPowerLogQuan(): The Quantile of the Power distribution must be computed correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogRand_ENABLED || setPowerLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i

        integer(IK) , parameter     :: NSIM = 1000_IK
        real(RKG)   , parameter     :: MEAN_REF = 0.5_RKG
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = 0.1_RKG
        real(RKG)                   :: alpha, logCDFNF, logMinX, logMaxX
        real(RKG)                   :: LogRand(NSIM), CDF(NSIM), mean, diff
        logical(LK)                 :: logMinXPresent

        assertion = .true._LK

        do i = 1_IK, 2_IK
            logMinX = getUnifRand(-LOG_HUGE, LOG_HUGE)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMinX = logMinX)
            call runTestsWith()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(logMinX)

            use pm_val2str, only: getStr
            use pm_sampleMean, only: getMean
            use pm_distNegExp, only: getNegExpRand
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(LOG_HUGE)
            real(RKG), intent(in), optional :: logMinX
            integer(IK) :: j

            logMinXPresent = present(logMinX)
            alpha = getUnifRand(0.5_RKG, +3._RKG)
            if (logMinXPresent) then
                logCDFNF = getPowerLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = getPowerLogCDFNF(alpha, logMaxX)
            end if

            LogRand = getUnifRand(getOption(logMaxX - 4._RKG, logMinX), logMaxX, s1 = size(LogRand, kind = IK))
            if (logMinXPresent) then
                do j = 1_IK, size(LogRand, kind = IK)
#if                 getPowerLogRand_ENABLED
                    LogRand(j) = getPowerLogRand(alpha, logMinX, logMaxX) ! Truncated Power distribution.
#elif               setPowerLogRand_ENABLED
                    call setPowerLogRand(LogRand(j), getNegExpRand(1._RKG), alpha, logMinX, logCDFNF)
#else
#error              "Unrecognized interface."
#endif
                    call setPowerLogCDF(CDF(j), LogRand(j), alpha, logMinX, logCDFNF)
                    call report(__LINE__, logMinX, LogRand(j))
                end do
            else
                do j = 1_IK, size(LogRand, kind = IK)
#if                 getPowerLogRand_ENABLED
                    LogRand(j) = getPowerLogRand(alpha, logMaxX) ! Truncated Power distribution.
#elif               setPowerLogRand_ENABLED
                    call setPowerLogRand(LogRand(j), getNegExpRand(1._RKG), alpha, logCDFNF)
#else
#error              "Unrecognized interface."
#endif
                    call setPowerLogCDF(CDF(j), LogRand(j), alpha, logCDFNF)
                    call report(__LINE__, logMinX, LogRand(j))
                end do
            end if
            mean = getMean(exp(CDF))
            call report(__LINE__, logMinX, mean = mean)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, logMinX, logRand, mean)
            integer, intent(in) :: line
            real(RKG), intent(in), optional :: logMinX, logRand, mean
            if (present(mean) .and. present(logRand)) then
                error stop "Internal error occurred. Both `logRand` and `mean` input arguments cannot be present simultaneously." ! LCOV_EXCL_LINE
            elseif (present(mean)) then
                diff = abs(mean - MEAN_REF)
                assertion = assertion .and. diff <= max(abs(MEAN_REF) * TOL, TOL)
            elseif (present(logRand)) then
                assertion = assertion .and. logical(logMaxX >= logRand, LK) .and. logical(logRand >= getOption(logMinX, -huge(logMinX)), LK)
            else
                error stop "Internal error occurred. Either `logRand` and `mean` input argument must be present." ! LCOV_EXCL_LINE
            end if
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(LogRand)", present(LogRand)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMinX)", present(logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, LogRand)", getOption(HUGE_RKG, LogRand)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMinX)", getOption(-HUGE_RKG, logMinX)
                write(test%disp%unit,"(*(g0,:,', '))") "logCDFNF", logCDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMaxX", logMaxX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha", alpha
                if (present(mean)) then
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "mean", mean
                write(test%disp%unit,"(*(g0,:,', '))") "MEAN_REF", MEAN_REF
                write(test%disp%unit,"(*(g0,:,', '))")
                end if
                ! LCOV_EXCL_STOP
            end if
            if (present(mean)) then
                call test%assert(assertion, SK_"The generated random value must be Power-distributed. Due to the statistical nature of the computations, this test can occasional fail, in which case, rerunning the test may fix the failure.", int(line, IK))
            elseif (present(logRand)) then
                call test%assert(assertion, SK_"The generated Power random value must be within the support of the distribution.", int(line, IK))
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif