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
!>  This include file contains procedure implementations of the tests of [pm_distPiwiPoweto](@ref pm_distPiwiPoweto).
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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getPiwiPowetoLogPDF_ENABLED || setPiwiPowetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, numComp

        real(RKC)   , parameter     :: TOL = sqrt(epsilon(0._RKC))
        real(RKC)   , parameter     :: LOG_HUGE = log(huge(0._RKC))
        real(RKC)   , allocatable   :: logLimX(:), alpha(:), logPDFNF(:)
        real(RKC)                   :: logx, output, output_ref, diff

        assertion = .true._LK

        do i = 1_IK, 20_IK

            call setUnifRand(numComp, 1_IK, 5_IK)
            !if (getUnifRand()) then
            !    logLimX = [getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp), LOG_HUGE]
            !else
            !    logLimX = getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp + 1)
            !end if
            !logLimX = getUnifRand(log(sqrt(epsilon(0._RKC))), sqrt(LOG_HUGE), numComp + 1)
            !logLimX = getUnifRand(log(sqrt(epsilon(0._RKC))), 10._RKC, numComp + 1)
            logLimX = getUnifRand(-5._RKC, 8._RKC, numComp + 1)
            call setSorted(logLimX)
            if (logLimX(size(logLimX)) == LOG_HUGE) then
                alpha = [getUnifRand(-2._RKC, 2._RKC), getUnifRand(-2._RKC, 2._RKC, numComp - 2_IK), getUnifRand(-8._RKC, -2._RKC)]
            else
                alpha = getUnifRand(-2._RKC, 2._RKC, numComp)
            end if
            if (getUnifRand() .and. size(alpha) > 1) alpha(getUnifRand(1, size(alpha) - 1)) = -1._RKC
            if (getUnifRand() .and. getUnifRand() .and. size(alpha) > 1) alpha(getUnifRand(1, size(alpha) - 1)) = 0._RKC
            if (getUnifRand()) then
                alpha(size(alpha)) = -1._RKC
            elseif (getUnifRand()) then
                alpha(size(alpha)) = 0._RKC
            end if
            logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getPiwiPowetoLogPDF_ENABLED

            logx = getUnifRand(logLimX(1), logLimX(size(logLimX)))
            output = getPiwiPowetoLogPDF(logx, alpha, logLimX)
            call setPiwiPowetoLogPDF(output_ref, logx, alpha, logLimX, logPDFNF)
            call report()
            call test%assert(assertion, SK_"@getPiwiPowetoLogPDF(): The PDF of the PiwiPoweto distribution must be computed correctly when `logPDFNF` is missing.", int(__LINE__, IK))

            logx = getUnifRand(logLimX(1), logLimX(size(logLimX)))
            output = getPiwiPowetoLogPDF(logx, alpha, logLimX, logPDFNF)
            call setPiwiPowetoLogPDF(output_ref, logx, alpha, logLimX, logPDFNF)
            call report()
            call test%assert(assertion, SK_"@getPiwiPowetoLogPDF(): The PDF of the PiwiPoweto distribution must be computed correctly when `logPDFNF` is present.", int(__LINE__, IK))

#elif       setPiwiPowetoLogPDF_ENABLED

            block

                use pm_distPiwiPoweto, only: getPiwiPowetoCDF
                use pm_arraySearch, only: getBin
                use pm_val2str, only: getStr
                character(255, SK) :: msg
                integer(IK) :: j

                !output_ref = 1._RKC
                msg = SK_" "
                !output = 0._RKC
                do j = 1, numComp - 1
                    output_ref = getPiwiPowetoCDF(logLimX(j + 1), alpha, logLimX) - getPiwiPowetoCDF(logLimX(j), alpha, logLimX)
                    assertion = assertion .and. .not. isFailedQuad(getPiwiPowetoPDF, exp(logLimX(j)), exp(logLimX(j + 1)), output, reltol = TOL, msg = msg)
                    call report()
                    call test%assert(assertion, SK_"@setPiwiPowetoLogPDF(): The integral of the PiwiPoweto distribution must be computed correctly.", int(__LINE__, IK))
                    !call test%assert(assertion, SK_"@setPiwiPowetoLogPDF(): The integration of the PiwiPoweto distribution must not fail. j, msg: "//getStr(j)//SK_", "//trim(msg), int(__LINE__, IK))
                    !output = output + logx
                end do

                !call report()
                !call test%assert(assertion, SK_"@setPiwiPowetoLogPDF(): The integral of the PiwiPoweto distribution must be computed correctly.", int(__LINE__, IK))

                logx = getUnifRand(logLimX(1), logLimX(size(logLimX)))
                call setPiwiPowetoLogPDF(output_ref, logx, alpha, logLimX, logPDFNF)
                call setPiwiPowetoLogPDF(output, logx, alpha, getBin(logLimX, logx), logPDFNF)
                call report()
                call test%assert(assertion, SK_"@setPiwiPowetoLogPDF(): The PDF of the PiwiPoweto distribution must be computed correctly when the value bin is specified.", int(__LINE__, IK))

            end block

#else
#error      "Unrecognized interface."
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setPiwiPowetoLogPDF_ENABLED
        function getPiwiPowetoPDF(x) result(pdf)
            real(RKC), intent(in) :: x
            real(RKC) :: pdf
            call setPiwiPowetoLogPDF(pdf, log(x), alpha, logLimX, logPDFNF)
            pdf = exp(pdf)
            !if (pdf > -huge(pdf)**0.9) then
            !    pdf = exp(pdf)
            !else
            !    pdf = 0._RKC
            !end if
        end function
#elif   !getPiwiPowetoLogPDF_ENABLED
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(numComp * abs(output_ref) * TOL, TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "numComp        ", numComp
                write(test%disp%unit,"(*(g0,:,', '))") "size(alpha)    ", size(alpha)
                write(test%disp%unit,"(*(g0,:,', '))") "size(logLimX)  ", size(logLimX)
                write(test%disp%unit,"(*(g0,:,', '))") "logPDFNF     ", logPDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logLimX        ", logLimX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha          ", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL            ", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output         ", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref     ", output_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPiwiPowetoCDF_ENABLED || setPiwiPowetoCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, numComp
        real(RKC)   , parameter     :: TOL = sqrt(epsilon(0._RKC))
        real(RKC)   , allocatable   :: logLimX(:), cumSumArea(:), alpha(:), logPDFNF(:)
        real(RKC)                   :: output, output_ref, diff, lb, ub

        assertion = .true._LK

        do i = 1_IK, 100_IK

            call setUnifRand(numComp, 1_IK, 5_IK)
            !if (getUnifRand()) then
            !    logLimX = [getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp), log(huge(0._RKC))]
            !else
            !    logLimX = getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp + 1)
            !end if
            logLimX = getUnifRand(-5._RKC, 5._RKC, numComp + 1)
            call setSorted(logLimX)
            if (numComp > 1_IK) then
                alpha = [getUnifRand(-1.01_RKC, 3._RKC), getUnifRand(-3._RKC, 3._RKC, numComp - 2_IK), getUnifRand(-8._RKC, -2._RKC)]
            else
                alpha = [getUnifRand(-8._RKC, -2._RKC)]
            end if
            if (getUnifRand() .and. size(alpha) > 1) alpha(getUnifRand(1, size(alpha) - 1)) = -1._RKC
            if (getUnifRand() .and. getUnifRand() .and. size(alpha) > 1) alpha(getUnifRand(1, size(alpha) - 1)) = 0._RKC
            if (size(logLimX) > size(alpha)) then
                if (getUnifRand()) then
                    alpha(size(alpha)) = -1._RKC
                elseif (getUnifRand()) then
                    alpha(size(alpha)) = 0._RKC
                end if
            end if
            if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea, mold = logLimX)
            logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea)

            ! ub = merge(logLimX(size(logLimX)), logLimX(size(logLimX)), size(logLimX,1,IK) > numComp)
            lb = getUnifRand(logLimX(1), logLimX(size(logLimX)))
            ub = exp(getUnifRand(lb, logLimX(size(logLimX))))
            lb = exp(lb)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getPiwiPowetoCDF_ENABLED

            block

                real(RKC) :: logx

                logx = merge(getUnifRand(logLimX(1), logLimX(size(logLimX))), getUnifRand(logLimX(1), logLimX(size(logLimX))), getUnifRand())
                output = getPiwiPowetoCDF(logx, alpha, logLimX)
                call setPiwiPowetoCDF(output_ref, logx, alpha, logLimX, logPDFNF, cumSumArea)
                call report()
                call test%assert(assertion, SK_"@getPiwiPowetoCDF(): The CDF of the PiwiPoweto distribution must be computed correctly when `logPDFNF` is missing.", int(__LINE__, IK))

                logx = merge(getUnifRand(logLimX(1), logLimX(size(logLimX))), getUnifRand(logLimX(1), logLimX(size(logLimX))), getUnifRand())
                output = getPiwiPowetoCDF(logx, alpha, logLimX, logPDFNF, cumSumArea)
                call setPiwiPowetoCDF(output_ref, logx, alpha, logLimX, logPDFNF, cumSumArea)
                call report()
                call test%assert(assertion, SK_"@getPiwiPowetoCDF(): The CDF of the PiwiPoweto distribution must be computed correctly when `logPDFNF` is present.", int(__LINE__, IK))

            end block

#elif       setPiwiPowetoCDF_ENABLED

            block

                use pm_arraySearch, only: getBin
                use pm_val2str, only: getStr
                real(RKC) :: upperCDF, lowerCDF
                character(255, SK) :: msg
                integer(IK) :: j

                ! Validate the CDF on the boundary points.

                do j = 1, numComp

                    output_ref = 0._RKC
                    if (j > 1_IK) output_ref = cumSumArea(j)

                    call setPiwiPowetoCDF(output, logLimX(j), alpha, logLimX, logPDFNF, cumSumArea)
                    call report()
                    call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The CDF at the component boundary values must equal the corresponding elements of `cumSumArea`.", int(__LINE__, IK))

                    call setPiwiPowetoCDF(output, logLimX(j), alpha, logLimX, logPDFNF, cumSumArea, bin = j)
                    call report()
                    call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The CDF at the component boundary values must equal the corresponding elements of `cumSumArea`.", int(__LINE__, IK))

                end do

                ! The CDF over the entire support must equal unity.

                output_ref = 1._RKC

                call setPiwiPowetoCDF(output, logLimX(size(logLimX)), alpha, logLimX, logPDFNF, cumSumArea)
                call report()
                call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The CDF over the entire support must equal unity.", int(__LINE__, IK))

                call setPiwiPowetoCDF(output, logLimX(size(logLimX)), alpha, logLimX, logPDFNF, cumSumArea, bin = numComp)
                call report()
                call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The CDF at the component boundary values must equal the corresponding elements of `cumSumArea` with `bin` present.", int(__LINE__, IK))

                msg = SK_" "

                ! Segments of CDF must be computed correctly.

                j = getBin(logLimX, log(lb))
                assertion = assertion .and. .not. isFailedQuad(getPiwiPowetoPDF, exp(logLimX(j)), lb, lowerCDF, reltol = TOL, msg = msg)
                call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The integration of the PiwiPoweto distribution must not fail. lb, ub, exp(logLimX), msg: "//getStr([lb, ub, exp(logLimX)])//SK_", "//trim(msg), int(__LINE__, IK))
                output_ref = cumSumArea(j) + lowerCDF

                j = getBin(logLimX, log(ub))
                assertion = assertion .and. .not. isFailedQuad(getPiwiPowetoPDF, exp(logLimX(j)), ub, upperCDF, reltol = TOL, msg = msg)
                call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The integration of the PiwiPoweto distribution must not fail. lb, ub, exp(logLimX), msg: "//getStr([lb, ub, exp(logLimX)])//SK_", "//trim(msg), int(__LINE__, IK))
                output_ref = cumSumArea(j) + upperCDF - output_ref

                call setPiwiPowetoCDF(lowerCDF, log(lb), alpha, logLimX, logPDFNF, cumSumArea)
                call setPiwiPowetoCDF(upperCDF, log(ub), alpha, logLimX, logPDFNF, cumSumArea)
                output = upperCDF - lowerCDF
                call report()
                call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The CDF between two arbitrary lower and upper limits must be computed correctly.", int(__LINE__, IK))

                call setPiwiPowetoCDF(lowerCDF, log(lb), alpha, logLimX, logPDFNF, cumSumArea, bin = getBin(logLimX, log(lb)))
                call setPiwiPowetoCDF(upperCDF, log(ub), alpha, logLimX, logPDFNF, cumSumArea, bin = getBin(logLimX, log(ub)))
                output = upperCDF - lowerCDF
                call report()
                call test%assert(assertion, SK_"@setPiwiPowetoCDF(): The CDF between two arbitrary lower and upper limits must be computed correctly with `bin` present.", int(__LINE__, IK))

            end block

#else
#error      "Unrecognized interface."
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setPiwiPowetoCDF_ENABLED
        function getPiwiPowetoPDF(x) result(pdf)
            use pm_distPiwiPoweto, only: setPiwiPowetoLogPDF
            real(RKC), intent(in) :: x
            real(RKC) :: pdf
            call setPiwiPowetoLogPDF(pdf, log(x), alpha, logLimX, logPDFNF)
            pdf = exp(pdf)
        end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(output - output_ref)
            assertion = assertion .and. diff <= max(numComp * abs(output_ref) * TOL, TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "numComp        ", numComp
                write(test%disp%unit,"(*(g0,:,', '))") "size(alpha)    ", size(alpha)
                write(test%disp%unit,"(*(g0,:,', '))") "size(logLimX)  ", size(logLimX)
                write(test%disp%unit,"(*(g0,:,', '))") "cumSumArea     ", cumSumArea
                write(test%disp%unit,"(*(g0,:,', '))") "logPDFNF     ", logPDFNF
                write(test%disp%unit,"(*(g0,:,', '))") "logLimX        ", logLimX
                write(test%disp%unit,"(*(g0,:,', '))") "alpha          ", alpha
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL           ", TOL
                write(test%disp%unit,"(*(g0,:,', '))") "output         ", output
                write(test%disp%unit,"(*(g0,:,', '))") "output_ref     ", output_ref
                write(test%disp%unit,"(*(g0,:,', '))") "lb, ub         ", lb, ub
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif