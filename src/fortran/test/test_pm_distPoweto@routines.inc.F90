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
!>  This include file contains procedure implementations of the tests of [pm_distPoweto](@ref pm_distPoweto).
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
#if     getPowetoLogPDF_ENABLED || setPowetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, numComp

        real(RKC)   , parameter     :: TOL = sqrt(epsilon(0._RKC))
        real(RKC)   , parameter     :: LOG_HUGE = log(huge(0._RKC))
        real(RKC)   , allocatable   :: LogLimX(:), Alpha(:), LogNormFac(:)
        real(RKC)                   :: logx, output, output_ref, diff

        assertion = .true._LK

        do i = 1_IK, 20_IK

            call setUnifRand(numComp, 1_IK, 5_IK)
            !if (getUnifRand()) then
            !    LogLimX = [getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp), LOG_HUGE]
            !else
            !    LogLimX = getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp + 1)
            !end if
            !LogLimX = getUnifRand(log(sqrt(epsilon(0._RKC))), sqrt(LOG_HUGE), numComp + 1)
            !LogLimX = getUnifRand(log(sqrt(epsilon(0._RKC))), 10._RKC, numComp + 1)
            LogLimX = getUnifRand(-5._RKC, 8._RKC, numComp + 1)
            call setSorted(LogLimX)
            if (LogLimX(size(LogLimX)) == LOG_HUGE) then
                Alpha = [getUnifRand(-2._RKC, 2._RKC), getUnifRand(-2._RKC, 2._RKC, numComp - 2_IK), getUnifRand(-8._RKC, -2._RKC)]
            else
                Alpha = getUnifRand(-2._RKC, 2._RKC, numComp)
            end if
            if (getUnifRand() .and. size(Alpha) > 1) Alpha(getUnifRand(1, size(Alpha) - 1)) = -1._RKC
            if (getUnifRand() .and. getUnifRand() .and. size(Alpha) > 1) Alpha(getUnifRand(1, size(Alpha) - 1)) = 0._RKC
            if (getUnifRand()) then
                Alpha(size(Alpha)) = -1._RKC
            elseif (getUnifRand()) then
                Alpha(size(Alpha)) = 0._RKC
            end if
            LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getPowetoLogPDF_ENABLED

            logx = getUnifRand(LogLimX(1), LogLimX(size(LogLimX)))
            output = getPowetoLogPDF(logx, Alpha, LogLimX)
            call setPowetoLogPDF(output_ref, logx, Alpha, LogLimX, LogNormFac)
            call report()
            call test%assert(assertion, SK_"@getPowetoLogPDF(): The PDF of the Poweto distribution must be computed correctly when `LogNormFac` is missing.", int(__LINE__, IK))

            logx = getUnifRand(LogLimX(1), LogLimX(size(LogLimX)))
            output = getPowetoLogPDF(logx, Alpha, LogLimX, LogNormFac)
            call setPowetoLogPDF(output_ref, logx, Alpha, LogLimX, LogNormFac)
            call report()
            call test%assert(assertion, SK_"@getPowetoLogPDF(): The PDF of the Poweto distribution must be computed correctly when `LogNormFac` is present.", int(__LINE__, IK))

#elif       setPowetoLogPDF_ENABLED

            block

                use pm_distPoweto, only: getPowetoCDF
                use pm_arraySearch, only: getBin
                use pm_val2str, only: getStr
                character(255, SK) :: msg
                integer(IK) :: j

                !output_ref = 1._RKC
                msg = SK_" "
                !output = 0._RKC
                do j = 1, numComp - 1
                    output_ref = getPowetoCDF(LogLimX(j + 1), Alpha, LogLimX) - getPowetoCDF(LogLimX(j), Alpha, LogLimX)
                    assertion = assertion .and. .not. isFailedQuad(getPowetoPDF, exp(LogLimX(j)), exp(LogLimX(j + 1)), output, reltol = TOL, msg = msg)
                    call report()
                    call test%assert(assertion, SK_"@setPowetoLogPDF(): The integral of the Poweto distribution must be computed correctly.", int(__LINE__, IK))
                    !call test%assert(assertion, SK_"@setPowetoLogPDF(): The integration of the Poweto distribution must not fail. j, msg: "//getStr(j)//SK_", "//trim(msg), int(__LINE__, IK))
                    !output = output + logx
                end do

                !call report()
                !call test%assert(assertion, SK_"@setPowetoLogPDF(): The integral of the Poweto distribution must be computed correctly.", int(__LINE__, IK))

                logx = getUnifRand(LogLimX(1), LogLimX(size(LogLimX)))
                call setPowetoLogPDF(output_ref, logx, Alpha, LogLimX, LogNormFac)
                call setPowetoLogPDF(output, logx, Alpha, getBin(LogLimX, logx), LogNormFac)
                call report()
                call test%assert(assertion, SK_"@setPowetoLogPDF(): The PDF of the Poweto distribution must be computed correctly when the value bin is specified.", int(__LINE__, IK))

            end block

#else
#error      "Unrecognized interface."
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setPowetoLogPDF_ENABLED
        function getPowetoPDF(x) result(pdf)
            real(RKC), intent(in) :: x
            real(RKC) :: pdf
            call setPowetoLogPDF(pdf, log(x), Alpha, LogLimX, LogNormFac)
            pdf = exp(pdf)
            !if (pdf > -huge(pdf)**0.9) then
            !    pdf = exp(pdf)
            !else
            !    pdf = 0._RKC
            !end if
        end function
#elif   !getPowetoLogPDF_ENABLED
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
                write(test%disp%unit,"(*(g0,:,', '))") "size(Alpha)    ", size(Alpha)
                write(test%disp%unit,"(*(g0,:,', '))") "size(LogLimX)  ", size(LogLimX)
                write(test%disp%unit,"(*(g0,:,', '))") "LogNormFac     ", LogNormFac
                write(test%disp%unit,"(*(g0,:,', '))") "LogLimX        ", LogLimX
                write(test%disp%unit,"(*(g0,:,', '))") "Alpha          ", Alpha
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
#elif   getPowetoCDF_ENABLED || setPowetoCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, numComp
        real(RKC)   , parameter     :: TOL = sqrt(epsilon(0._RKC))
        real(RKC)   , allocatable   :: LogLimX(:), CumSumArea(:), Alpha(:), LogNormFac(:)
        real(RKC)                   :: output, output_ref, diff, lb, ub

        assertion = .true._LK

        do i = 1_IK, 100_IK

            call setUnifRand(numComp, 1_IK, 5_IK)
            !if (getUnifRand()) then
            !    LogLimX = [getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp), log(huge(0._RKC))]
            !else
            !    LogLimX = getUnifRand(sqrt(epsilon(0._RKC)), 10._RKC, numComp + 1)
            !end if
            LogLimX = getUnifRand(-5._RKC, 5._RKC, numComp + 1)
            call setSorted(LogLimX)
            if (numComp > 1_IK) then
                Alpha = [getUnifRand(-1.01_RKC, 3._RKC), getUnifRand(-3._RKC, 3._RKC, numComp - 2_IK), getUnifRand(-8._RKC, -2._RKC)]
            else
                Alpha = [getUnifRand(-8._RKC, -2._RKC)]
            end if
            if (getUnifRand() .and. size(Alpha) > 1) Alpha(getUnifRand(1, size(Alpha) - 1)) = -1._RKC
            if (getUnifRand() .and. getUnifRand() .and. size(Alpha) > 1) Alpha(getUnifRand(1, size(Alpha) - 1)) = 0._RKC
            if (size(LogLimX) > size(Alpha)) then
                if (getUnifRand()) then
                    Alpha(size(Alpha)) = -1._RKC
                elseif (getUnifRand()) then
                    Alpha(size(Alpha)) = 0._RKC
                end if
            end if
            if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea, mold = LogLimX)
            LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea)

            ! ub = merge(LogLimX(size(LogLimX)), LogLimX(size(LogLimX)), size(LogLimX,1,IK) > numComp)
            lb = getUnifRand(LogLimX(1), LogLimX(size(LogLimX)))
            ub = exp(getUnifRand(lb, LogLimX(size(LogLimX))))
            lb = exp(lb)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getPowetoCDF_ENABLED

            block

                real(RKC) :: logx

                logx = merge(getUnifRand(LogLimX(1), LogLimX(size(LogLimX))), getUnifRand(LogLimX(1), LogLimX(size(LogLimX))), getUnifRand())
                output = getPowetoCDF(logx, Alpha, LogLimX)
                call setPowetoCDF(output_ref, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
                call report()
                call test%assert(assertion, SK_"@getPowetoCDF(): The CDF of the Poweto distribution must be computed correctly when `LogNormFac` is missing.", int(__LINE__, IK))

                logx = merge(getUnifRand(LogLimX(1), LogLimX(size(LogLimX))), getUnifRand(LogLimX(1), LogLimX(size(LogLimX))), getUnifRand())
                output = getPowetoCDF(logx, Alpha, LogLimX, LogNormFac, CumSumArea)
                call setPowetoCDF(output_ref, logx, Alpha, LogLimX, LogNormFac, CumSumArea)
                call report()
                call test%assert(assertion, SK_"@getPowetoCDF(): The CDF of the Poweto distribution must be computed correctly when `LogNormFac` is present.", int(__LINE__, IK))

            end block

#elif       setPowetoCDF_ENABLED

            block

                use pm_arraySearch, only: getBin
                use pm_val2str, only: getStr
                real(RKC) :: upperCDF, lowerCDF
                character(255, SK) :: msg
                integer(IK) :: j

                ! Validate the CDF on the boundary points.

                do j = 1, numComp

                    output_ref = 0._RKC
                    if (j > 1_IK) output_ref = CumSumArea(j)

                    call setPowetoCDF(output, LogLimX(j), Alpha, LogLimX, LogNormFac, CumSumArea)
                    call report()
                    call test%assert(assertion, SK_"@setPowetoCDF(): The CDF at the component boundary values must equal the corresponding elements of `CumSumArea`.", int(__LINE__, IK))

                    call setPowetoCDF(output, LogLimX(j), Alpha, LogLimX, LogNormFac, CumSumArea, bin = j)
                    call report()
                    call test%assert(assertion, SK_"@setPowetoCDF(): The CDF at the component boundary values must equal the corresponding elements of `CumSumArea`.", int(__LINE__, IK))

                end do

                ! The CDF over the entire support must equal unity.

                output_ref = 1._RKC

                call setPowetoCDF(output, LogLimX(size(LogLimX)), Alpha, LogLimX, LogNormFac, CumSumArea)
                call report()
                call test%assert(assertion, SK_"@setPowetoCDF(): The CDF over the entire support must equal unity.", int(__LINE__, IK))

                call setPowetoCDF(output, LogLimX(size(LogLimX)), Alpha, LogLimX, LogNormFac, CumSumArea, bin = numComp)
                call report()
                call test%assert(assertion, SK_"@setPowetoCDF(): The CDF at the component boundary values must equal the corresponding elements of `CumSumArea` with `bin` present.", int(__LINE__, IK))

                msg = SK_" "

                ! Segments of CDF must be computed correctly.

                j = getBin(LogLimX, log(lb))
                assertion = assertion .and. .not. isFailedQuad(getPowetoPDF, exp(LogLimX(j)), lb, lowerCDF, reltol = TOL, msg = msg)
                call test%assert(assertion, SK_"@setPowetoCDF(): The integration of the Poweto distribution must not fail. lb, ub, exp(LogLimX), msg: "//getStr([lb, ub, exp(LogLimX)])//SK_", "//trim(msg), int(__LINE__, IK))
                output_ref = CumSumArea(j) + lowerCDF

                j = getBin(LogLimX, log(ub))
                assertion = assertion .and. .not. isFailedQuad(getPowetoPDF, exp(LogLimX(j)), ub, upperCDF, reltol = TOL, msg = msg)
                call test%assert(assertion, SK_"@setPowetoCDF(): The integration of the Poweto distribution must not fail. lb, ub, exp(LogLimX), msg: "//getStr([lb, ub, exp(LogLimX)])//SK_", "//trim(msg), int(__LINE__, IK))
                output_ref = CumSumArea(j) + upperCDF - output_ref

                call setPowetoCDF(lowerCDF, log(lb), Alpha, LogLimX, LogNormFac, CumSumArea)
                call setPowetoCDF(upperCDF, log(ub), Alpha, LogLimX, LogNormFac, CumSumArea)
                output = upperCDF - lowerCDF
                call report()
                call test%assert(assertion, SK_"@setPowetoCDF(): The CDF between two arbitrary lower and upper limits must be computed correctly.", int(__LINE__, IK))

                call setPowetoCDF(lowerCDF, log(lb), Alpha, LogLimX, LogNormFac, CumSumArea, bin = getBin(LogLimX, log(lb)))
                call setPowetoCDF(upperCDF, log(ub), Alpha, LogLimX, LogNormFac, CumSumArea, bin = getBin(LogLimX, log(ub)))
                output = upperCDF - lowerCDF
                call report()
                call test%assert(assertion, SK_"@setPowetoCDF(): The CDF between two arbitrary lower and upper limits must be computed correctly with `bin` present.", int(__LINE__, IK))

            end block

#else
#error      "Unrecognized interface."
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setPowetoCDF_ENABLED
        function getPowetoPDF(x) result(pdf)
            use pm_distPoweto, only: setPowetoLogPDF
            real(RKC), intent(in) :: x
            real(RKC) :: pdf
            call setPowetoLogPDF(pdf, log(x), Alpha, LogLimX, LogNormFac)
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
                write(test%disp%unit,"(*(g0,:,', '))") "size(Alpha)    ", size(Alpha)
                write(test%disp%unit,"(*(g0,:,', '))") "size(LogLimX)  ", size(LogLimX)
                write(test%disp%unit,"(*(g0,:,', '))") "CumSumArea     ", CumSumArea
                write(test%disp%unit,"(*(g0,:,', '))") "LogNormFac     ", LogNormFac
                write(test%disp%unit,"(*(g0,:,', '))") "LogLimX        ", LogLimX
                write(test%disp%unit,"(*(g0,:,', '))") "Alpha          ", Alpha
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