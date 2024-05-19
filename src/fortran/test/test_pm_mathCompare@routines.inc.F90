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
!>  This include file contains procedure implementations of the tests of [pm_mathCompare](@ref pm_mathCompare).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_distUnif, only: getUnifRand, setUnifRand
        use pm_except, only: setInfNeg, setInfPos, isInf
        use pm_except, only: setNAN, isNAN
#if     isClose_CK_ENABLED
        complex(CKG)                    :: x, y
        complex(CKG)    , parameter     :: ZERO = (0._CKG, 0._CKG)
        complex(CKG)    , parameter     :: low = -cmplx(sqrt(sqrt(huge(0._CKG))), sqrt(sqrt(huge(0._CKG))), CKG)
        complex(CKG)    , parameter     :: upp = +cmplx(sqrt(sqrt(huge(0._CKG))), sqrt(sqrt(huge(0._CKG))), CKG)
#elif   isClose_RK_ENABLED
        real(RKG)                       :: x, y
        real(RKG)       , parameter     :: low = -sqrt(sqrt(huge(low)))
        real(RKG)       , parameter     :: upp = +sqrt(sqrt(huge(upp)))
        real(RKG)       , parameter     :: ZERO = 0._RKG
#else
#error  "Unrecognized interface."
#endif
        integer         , parameter     :: TKG = kind(x)
        real(TKG)       , parameter     :: EPS = epsilon(1._TKG) * 100, TIN = tiny(0._TKG)
        logical(LK)                     :: close, close_def, isRelTolDef, isAbsTolDef
        real(TKG)                       :: reltol, abstol
        class(*)        , allocatable   :: method
        integer(IK)                     :: i, choice

        reltol = EPS
        abstol = TIN
        assertion = .true._LK
        do i = 1, 1000

            choice = getUnifRand(0_IK, 9_IK)
            if (choice == 0_IK) then
                call setNAN(x)
            elseif (choice == 1_IK) then
                call setInfNeg(x)
            elseif (choice == 2_IK) then
                call setInfPos(x)
            elseif (choice == 3_IK) then
                x = ZERO
            else
                call setUnifRand(x, low, upp)
            end if

            isAbsTolDef = getUnifRand()
            isRelTolDef = getUnifRand()
            if (.not. (isInf(x) .or. isNAN(x))) reltol = merge(EPS, getUnifRand(0.1_TKG * abs(x), 10._TKG * abs(x)), isRelTolDef)
            if (.not. (isInf(x) .or. isNAN(x))) abstol = merge(TIN, getUnifRand(0.2_TKG * abs(x),  5._TKG * abs(x)), isAbsTolDef)

            choice = getUnifRand(0_IK, 9_IK)
            if (choice == 0_IK) then
                call setNAN(y)
            elseif (choice == 1_IK) then
                call setInfNeg(y)
            elseif (choice == 2_IK) then
                call setInfPos(y)
            elseif (choice == 3_IK) then
                y = ZERO
            elseif (.not. (isInf(x) .or. isNAN(x))) then
                call setUnifRand(y, x - 3 * reltol, x + 3 * reltol)
            else
                call setUnifRand(y, ZERO - 3 * reltol, ZERO + 3 * reltol)
            end if

            if (getUnifRand()) then ! method
                choice = getUnifRand(1_IK, 4_IK)
                if (choice == 1_IK) then
                    method = reference_type()
                elseif (choice == 2_IK) then
                    method = strong_type()
                elseif (choice == 3_IK) then
                    method = weak_type()
                elseif (choice == 4_IK) then
                    method = mean_type()
                end if
                if (getUnifRand()) then ! reltol
                    if (getUnifRand()) then ! abstol
                        call runTestsWith(method, reltol, abstol)
                    else ! no abstol
                        call runTestsWith(method, reltol)
                    end if
                else ! no reltol
                    if (getUnifRand()) then ! abstol
                        call runTestsWith(method, abstol = abstol)
                    else ! no abstol
                        call runTestsWith(method)
                    end if
                end if
            else ! no method
                if (getUnifRand()) then ! reltol
                    if (getUnifRand()) then ! abstol
                        call runTestsWith(reltol = reltol, abstol = abstol)
                    else ! no abstol
                        call runTestsWith(reltol = reltol)
                    end if
                else ! no reltol
                    if (getUnifRand()) then ! abstol
                        call runTestsWith(abstol = abstol)
                    else ! no abstol
                        call runTestsWith()
                    end if
                end if
            endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(method, reltol, abstol)
            class(*)    , intent(in), optional  :: method
            real(TKG)   , intent(in), optional  :: reltol, abstol
            if (present(method)) then
                close_def = isClose_def(method, reltol, abstol)
                select type (method)
                    type is (reference_type)
                        close = isClose(x, y, method, reltol, abstol)
                    type is (strong_type)
                        close = isClose(x, y, method, reltol, abstol)
                    type is (weak_type)
                        close = isClose(x, y, method, reltol, abstol)
                    type is (mean_type)
                        close = isClose(x, y, method, reltol, abstol)
                end select
            else
                close_def = isClose_def(reltol = reltol, abstol = abstol)
                close = isClose(x, y, reltol, abstol)
            endif
            call report(__LINE__, method, reltol, abstol)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isClose_def(method, reltol, abstol) result(close)
            use pm_except, only: isInf, isInfNeg, isInfPos
            use pm_except, only: isNAN
            class(*)    , intent(in), optional  :: method
            real(TKG)   , intent(in), optional  :: reltol, abstol
            logical(LK) :: close

            real(TKG) :: reltol_def, abstol_def
            real(TKG) :: absDiff

            if (isNAN(x) .or. isNAN(y)) then
                close = .false._LK
                return
            end if

            if (.not. (isInf(x) .or. isInf(y))) then
                close = logical(x == y, LK)
                if (close) return
                if (present(reltol)) then
                    reltol_def = reltol
                else
                    reltol_def = epsilon(0._TKG)
                end if
                if (present(abstol)) then
                    abstol_def = abstol
                else
                    abstol_def = tiny(0._TKG)
                end if
                absDiff = abs(y - x)
                if (present(method)) then
                    select type (method)
                        type is (reference_type)
                            close = logical(absDiff <= abs(reltol_def * x) .or. absDiff <= abstol_def, LK)
                        type is (strong_type)
                            close = logical((absDiff <= abs(reltol_def * x) .and. absDiff <= abs(reltol_def * y)) .or. absDiff <= abstol_def, LK)
                        type is (weak_type)
                            close = logical(absDiff <= abs(reltol_def * x) .or. absDiff <= abs(reltol_def * y) .or. absDiff <= abstol_def, LK)
                        type is (mean_type)
                            close = logical(absDiff <= abs(reltol_def * 0.5_TKG * (x + y)) .or. absDiff <= abstol_def, LK)
                    end select
                else ! assume `weak_type`.
                    close = logical(absDiff <= abs(reltol_def * x) .or. absDiff <= abs(reltol_def * y) .or. absDiff <= abstol_def, LK)
                end if
                return
            end if
            close = (isInfNeg(x) .and. isInfNeg(y)) .or. (isInfPos(x) .and. isInfPos(y))
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, method, reltol, abstol)
            use pm_io, only: display_type
            use pm_val2str, only: getStr
            use pm_option, only: getOption
            integer     , intent(in) :: line
            class(*)    , intent(in), optional :: method
            real(TKG)   , intent(in), optional :: reltol, abstol
            type(display_type)  :: disp
            assertion = assertion .and. (close .eqv. close_def)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip()
                call disp%show("[x, y]")
                call disp%show( [x, y] )
                call disp%show("[present(reltol), present(abstol)]")
                call disp%show( [present(reltol), present(abstol)] )
                call disp%show("[getOption(EPS, reltol), getOption(TIN, abstol)]")
                call disp%show( [getOption(EPS, reltol), getOption(TIN, abstol)] )
                call disp%show("present(method)")
                call disp%show( present(method) )
                if (present(method)) then
                    select type (method)
                        type is (reference_type)
                            call disp%show("method: reference_type()")
                        type is (strong_type)
                            call disp%show("method: strong_type()")
                        type is (weak_type)
                            call disp%show("method: weak_type()")
                        type is (mean_type)
                            call disp%show("method: mean_type()")
                    end select
                end if
                call disp%show("[close, close_def]")
                call disp%show( [close, close_def] )
                call disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The proximity of the two numbers must be computed correctly.", int(line, IK))
        end subroutine
