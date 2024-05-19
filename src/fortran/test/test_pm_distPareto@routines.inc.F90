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
!>  This include file contains procedure implementations of the tests of [pm_distPareto](@ref pm_distPareto).
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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getParetoLogPDF_ENABLED || setParetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = sqrt(epsilon(0._RKG))
        real(RKG)                   :: logx, alpha, logPDFNF, logMinX, logMaxX
        real(RKG)                   :: output, output_ref, diff

        assertion = .true._LK

        do i = 1_IK, 20_IK
            logMinX = getUnifRand(-5._RKG, +4._RKG)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMaxX)
            call runTestsWith()
        end do

    contains

        subroutine runTestsWith(logMaxX)
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(log(huge(0._RKG)))
            real(RKG), intent(in), optional :: logMaxX
            alpha = -getUnifRand(0.5_RKG, +3._RKG)
            if (present(logMaxX)) then
                logPDFNF = getParetoLogPDFNF(alpha, logMinX, logMaxX)
            else
                logPDFNF = getParetoLogPDFNF(alpha, logMinX)
            end if
#if         getParetoLogPDF_ENABLED
            logx = getUnifRand(logMinX, getOption(SQRT_LOG_HUGE, logMaxX))
            call setParetoLogPDF(output_ref, logx, alpha, logPDFNF)
            if (present(logMaxX)) then
                output = getParetoLogPDF(logx, alpha, logMinX, logMaxX)
            else
                output = getParetoLogPDF(logx, alpha, logMinX)
            end if
            call compare(__LINE__, logMaxX, logx)
#elif       setParetoLogPDF_ENABLED
            output_ref = 1._RKG
            block
                character(255, SK) :: msg
                msg = SK_" "
                logx = huge(0._RKG)
                if (present(logMaxX)) logx = exp(logMaxX)
                assertion = assertion .and. .not. isFailedQuad(getParetoPDF, exp(logMinX), logx, output, reltol = TOL, msg = msg)
                if (.not. present(logMaxX)) logx = LOG_HUGE
                call report(logMaxX, logx)
                call test%assert(assertion, SK_"@setParetoLogPDF(): The integral of the Pareto distribution must be computed without failure. msg: "//trim(msg), int(__LINE__, IK))
            end block
            call compare()
            call report(logMaxX)
            call test%assert(assertion, SK_"@setParetoLogPDF(): The PDF of the Pareto distribution must be computed correctly.", int(__LINE__, IK))
#else
#error      "Unrecognized interface."
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getParetoLogPDF_ENABLED
        subroutine compare(line, logMaxX, logx)
            integer, intent(in) :: line
            real(RKG), intent(in), optional :: logMaxX, logx
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
            call report(logMaxX, logx)
            call test%assert(assertion, SK_"@getParetoLogPDF(): The PDF of the Pareto distribution must be computed correctly.", int(line, IK))
        end subroutine
#elif   setParetoLogPDF_ENABLED
        subroutine compare()
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
        end subroutine
        function getParetoPDF(x) result(pdf)
            real(RKG), intent(in) :: x
            real(RKG) :: pdf
            call setParetoLogPDF(pdf, log(x), alpha, logPDFNF)
            pdf = exp(pdf)
        end function
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(logMaxX, logx)
            real(RKG), intent(in), optional :: logMaxX, logx
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(logx)", present(logx)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMaxX)", present(logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, logx)", getOption(HUGE_RKG, logx)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMaxX)", getOption(-HUGE_RKG, logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "logPDFNF", logPDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMinX", logMinX
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
#elif   getParetoLogCDF_ENABLED || setParetoLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i

        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = epsilon(0._RKG) * 10000
        real(RKG)                   :: alpha, logCDFNF, logMinX, logMaxX
        real(RKG)                   :: output, output_ref, diff
        logical(LK)                 :: logMaxXPresent

        assertion = .true._LK

        do i = 1_IK, 300_IK
            logMinX = getUnifRand(-5._RKG, +4._RKG)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMaxX = logMaxX)
            call runTestsWith()
        end do

    contains

        subroutine runTestsWith(logMaxX)

            use pm_val2str, only: getStr
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(LOG_HUGE)
            real(RKG), intent(in), optional :: logMaxX
            real(RKG) :: logx
            logMaxXPresent = present(logMaxX)
            alpha = -getUnifRand(sqrt(epsilon(0._RKG)), +3._RKG)
            if (logMaxXPresent) then
                logCDFNF = getParetoLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = 0._RKG
            end if

#if         getParetoLogCDF_ENABLED

            logx = getUnifRand(logMinX, getOption(logMinX + 4._RKG, logMaxX))
            if (logMaxXPresent) then
                call setParetoLogCDF(output_ref, logx, alpha, logMinX, logCDFNF)
                output = getParetoLogCDF(logx, alpha, logMinX, logMaxX)
            else
                call setParetoLogCDF(output_ref, logx, alpha, logMinX)
                output = getParetoLogCDF(logx, alpha, logMinX)
            end if
            call compare(__LINE__, logMaxX, logx)

#elif       setParetoLogCDF_ENABLED

            output_ref = -LOG_HUGE
            logx = logMinX
            if (logMaxXPresent) then
                call setParetoLogCDF(output, logx, alpha, logMinX, logCDFNF)
            else
                call setParetoLogCDF(output, logx, alpha, logMinX)
            end if
            call report(logMaxX, logx)
            call test%assert(assertion, SK_"@setParetoLogCDF(): The CDF of the Pareto distribution at the lower limit of the support must be zero.", int(__LINE__, IK))

            output_ref = 0._RKG
            logx = getOption(+LOG_HUGE, logMaxX)
            if (logMaxXPresent) then
                call setParetoLogCDF(output, logx, alpha, logMinX, logCDFNF)
            else
                call setParetoLogCDF(output, logx, alpha, logMinX)
            end if
            call report(logMaxX, logx)
            call test%assert(assertion, SK_"@setParetoLogCDF(): The CDF of the Pareto distribution at the upper limit of the support must be unity.", int(__LINE__, IK))

            ! Compare the CDF difference at random points with the corresponding integral of the PDF.

            block
                character(255, SK) :: msg
                msg = SK_" "
                logx = getUnifRand(logMinX, getOption(logMinX + 4._RKG, logMaxX))
                assertion = assertion .and. .not. isFailedQuad(getParetoPDF, exp(logMinX), exp(logx), output_ref, reltol = TOL, msg = msg)
                call report(logMaxX, logx)
                call test%assert(assertion, SK_"@setParetoLogCDF(): The integral of the Pareto distribution must be computed without failure. msg: "//trim(msg), int(__LINE__, IK))
                if (logMaxXPresent) then
                    call setParetoLogCDF(output, logx, alpha, logMinX, logCDFNF)
                else
                    call setParetoLogCDF(output, logx, alpha, logMinX)
                end if
                output = exp(output)
                call compare()
                call report(logMaxX, logx)
                call test%assert(assertion, SK_"@setParetoLogCDF(): The CDF of the Pareto distribution must be computed correctly.", int(__LINE__, IK))
            end block

#else
#error      "Unrecognized interface."
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getParetoLogCDF_ENABLED
        subroutine compare(line, logMaxX, logx)
            integer :: line
            real(RKG), intent(in), optional :: logMaxX, logx
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL) * 10
            call report(logMaxX, logx)
            call test%assert(assertion, SK_"@getParetoLogCDF(): The CDF of the Pareto distribution must be computed correctly.", int(line, IK))
        end subroutine
#elif   setParetoLogCDF_ENABLED
        subroutine compare()
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL) * 10
        end subroutine
        function getParetoPDF(x) result(pdf)
            use pm_distPareto, only: getParetoLogPDF
            real(RKG), intent(in) :: x
            real(RKG) :: pdf
            if (logMaxXPresent) then
                pdf = exp(getParetoLogPDF(log(x), alpha, logMinX, logMaxX))
            else
                pdf = exp(getParetoLogPDF(log(x), alpha, logMinX))
            end if
        end function
#else
#error  "Unrecognized interface."
#endif
        subroutine report(logMaxX, logx)
            real(RKG), intent(in), optional :: logMaxX, logx
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(logx)", present(logx)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMaxX)", present(logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, logx)", getOption(HUGE_RKG, logx)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMaxX)", getOption(-HUGE_RKG, logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "logCDFNF", logCDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMinX", logMinX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref", output_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogQuan_ENABLED || setParetoLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = sqrt(epsilon(0._RKG)) * 1000
        real(RKG)                   :: alpha, logCDFNF, logMinX, logMaxX
        real(RKG)                   :: logx, output, output_ref, diff
        logical(LK)                 :: logMaxXPresent

        assertion = .true._LK

        do i = 1_IK, 500_IK
            logMinX = getUnifRand(-LOG_HUGE, LOG_HUGE)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMaxX = logMaxX)
            call runTestsWith()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(logMaxX)

            use pm_val2str, only: getStr
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(LOG_HUGE)
            real(RKG), intent(in), optional :: logMaxX

            logMaxXPresent = present(logMaxX)
            alpha = -getUnifRand(0.5_RKG, +3._RKG)
            if (logMaxXPresent) then
                logCDFNF = getParetoLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = 0._RKG ! getParetoLogCDFNF(alpha, logMinX)
            end if


            logx = getUnifRand(logMinX, getOption(logMinX + 4._RKG, logMaxX))
            if (logMaxXPresent) then
                call setParetoLogCDF(output_ref, logx, alpha, logMinX, logCDFNF)
                if (output_ref > 0._RKG) output_ref = 0._RKG ! ensure cdf cannot be larger than 1 due to roundoff error, otherwise algorithm checks will catch it.
                !write(*,*) output_ref, logx, alpha, logMinX
#if             getParetoLogQuan_ENABLED
                output = getParetoLogQuan(output_ref, alpha, logMinX, logMaxX)
#elif           setParetoLogQuan_ENABLED
                call setParetoLogQuan(output, output_ref, alpha, logMinX, logCDFNF)
#else
#error          "Unrecognized interface."
#endif
            else
                call setParetoLogCDF(output_ref, logx, alpha, logMinX)
                !write(*,*) output_ref, logx, alpha, logMinX
#if             getParetoLogQuan_ENABLED
                output = getParetoLogQuan(output_ref, alpha, logMinX)
#elif           setParetoLogQuan_ENABLED
                call setParetoLogQuan(output, output_ref, alpha, logMinX)
#endif
            end if
            output_ref = logx
            call report(__LINE__, logMaxX, logx)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, logMaxX, logx)
            integer, intent(in) :: line
            real(RKG), intent(in), optional :: logMaxX, logx
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(abs(output_ref) * TOL, TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(logx)", present(logx)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMaxX)", present(logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, logx)", getOption(HUGE_RKG, logx)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMaxX)", getOption(-HUGE_RKG, logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "logCDFNF", logCDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMinX", logMinX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref", output_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@getParetoLogQuan(): The Quantile of the Pareto distribution must be computed correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogRand_ENABLED || setParetoLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        integer(IK) , parameter     :: NSIM = 1000_IK
        real(RKG)   , parameter     :: MEAN_REF = 0.5_RKG
        real(RKG)   , parameter     :: HUGE_RKG = huge(0._RKG)
        real(RKG)   , parameter     :: LOG_HUGE = log(HUGE_RKG)
        real(RKG)   , parameter     :: TOL = 0.1_RKG
        real(RKG)                   :: alpha, logCDFNF, logMinX, logMaxX
        real(RKG)                   :: LogRand(NSIM), CDF(NSIM), mean, diff
        logical(LK)                 :: logMaxXPresent

        assertion = .true._LK

        do i = 1_IK, 2_IK
            logMinX = getUnifRand(-LOG_HUGE, LOG_HUGE)
            logMaxX = getUnifRand(logMinX + 1._RKG, logMinX + 4._RKG)
            call runTestsWith(logMaxX = logMaxX)
            call runTestsWith()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(logMaxX)

            use pm_val2str, only: getStr
            use pm_sampleMean, only: getMean
            use pm_distNegExp, only: getNegExpRand
            real(RKG), parameter :: SQRT_LOG_HUGE = sqrt(LOG_HUGE)
            real(RKG), intent(in), optional :: logMaxX
            integer(IK) :: j

            logMaxXPresent = present(logMaxX)
            alpha = -getUnifRand(0.5_RKG, +3._RKG)
            if (logMaxXPresent) then
                logCDFNF = getParetoLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = 0._RKG ! getParetoLogCDFNF(alpha, logMinX)
            end if

            LogRand = getUnifRand(logMinX, getOption(logMinX + 4._RKG, logMaxX), s1 = size(LogRand, kind = IK))
            if (logMaxXPresent) then
                do j = 1_IK, size(LogRand, kind = IK)
#if                 getParetoLogRand_ENABLED
                    LogRand(j) = getParetoLogRand(alpha, logMinX, logMaxX) ! Truncated Pareto distribution.
#elif               setParetoLogRand_ENABLED
                    call setParetoLogRand(LogRand(j), getNegExpRand(1._RKG), alpha, logMinX, logCDFNF)
#else
#error              "Unrecognized interface."
#endif
                    call setParetoLogCDF(CDF(j), LogRand(j), alpha, logMinX, logCDFNF)
                    call report(__LINE__, logMaxX, LogRand(j))
                end do
            else
                do j = 1_IK, size(LogRand, kind = IK)
#if                 getParetoLogRand_ENABLED
                    LogRand(j) = getParetoLogRand(alpha, logMinX) ! Truncated Pareto distribution.
#elif               setParetoLogRand_ENABLED
                    call setParetoLogRand(LogRand(j), getNegExpRand(1._RKG), alpha, logMinX)
#else
#error              "Unrecognized interface."
#endif
                    call setParetoLogCDF(CDF(j), LogRand(j), alpha, logMinX)
                    call report(__LINE__, logMaxX, LogRand(j))
                end do
            end if
            mean = getMean(exp(CDF))
            call report(__LINE__, logMaxX, mean = mean)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, logMaxX, logRand, mean)
            integer, intent(in) :: line
            real(RKG), intent(in), optional :: logMaxX, logRand, mean
            if (present(mean) .and. present(logRand)) then
                error stop "Internal error occurred. Both `logRand` and `mean` input arguments cannot be present simultaneously." ! LCOV_EXCL_LINE
            elseif (present(mean)) then
                diff = abs(mean - MEAN_REF)
                assertion = assertion .and. diff <= max(abs(MEAN_REF) * TOL, TOL)
            elseif (present(logRand)) then
                assertion = assertion .and. logical(logMinX <= logRand, LK) .and. logical(logRand <= getOption(logMaxX, huge(logMaxX)), LK)
            else
                error stop "Internal error occurred. Either `logRand` and `mean` input argument must be present." ! LCOV_EXCL_LINE
            end if
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "present(LogRand)", present(LogRand)
                write(test%disp%unit,"(*(g0,:,', '))") "present(logMaxX)", present(logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(HUGE_RKG, LogRand)", getOption(HUGE_RKG, LogRand)
                write(test%disp%unit,"(*(g0,:,', '))") "getOption(-HUGE_RKG, logMaxX)", getOption(-HUGE_RKG, logMaxX)
                write(test%disp%unit,"(*(g0,:,', '))") "logCDFNF", logCDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logMinX", logMinX
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
                call test%assert(assertion, SK_"The generated random value must be Pareto-distributed. Due to the statistical nature of the computations, this test can occasional fail, in which case, rerunning the test may fix the failure.", int(line, IK))
            elseif (present(logRand)) then
                call test%assert(assertion, SK_"The generated Pareto random value must be within the support of the distribution.", int(line, IK))
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif