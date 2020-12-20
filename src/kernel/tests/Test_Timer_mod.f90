!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [Timer_mod](@ref timer_mod).
!>  \author Amir Shahmoradi

module Test_Timer_mod

    use Timer_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Timer

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Timer()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_Timer_type_1, "test_Timer_type_1")
        call Test%run(test_getTimeSinceStart_1, "test_getTimeSinceStart_1")
        call Test%run(test_getTimeSinceLastCall_1, "test_getTimeSinceLastCall_1")
        call Test%finalize()

    end subroutine test_Timer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_Timer_type_1() result(assertion)

        use System_mod, only: sleep
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), parameter :: seconds = 0.05_RK
        logical             :: assertion
        type(Timer_type)    :: Timer
        type(Err_type)      :: Err

        Timer = Timer_type(Err)
        assertion = .not. Err%occurred; if (.not. assertion) return
        call sleep(seconds=seconds,Err=Err)
        assertion = .not. Err%occurred; if (.not. assertion) return
        call Timer%toc()
        assertion = assertion .and. Timer%Time%total > 0._RK
        assertion = assertion .and. Timer%Time%delta > 0._RK
        assertion = assertion .and. Timer%Time%start < Timer%Time%stop
        assertion = assertion .and. Timer%period > 0._RK

        assertion = assertion .and. Timer%Count%total > 0_IK
        assertion = assertion .and. Timer%Count%delta > 0_IK
        assertion = assertion .and. Timer%Count%start < Timer%Count%stop
        assertion = assertion .and. Timer%Count%rate > 0_IK
        assertion = assertion .and. Timer%Count%max > 0_IK

        assertion = assertion .and. Timer%Time%delta > 0.9_RK * seconds

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%start: ", Timer%Count%start
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%rate : ", Timer%Count%rate
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%max  : ", Timer%Count%max
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%stop : ", Timer%Count%stop
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%total: ", Timer%Count%total
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%delta: ", Timer%Count%delta
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%start : ", Timer%Time%start
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%stop  : ", Timer%Time%stop
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%total : ", Timer%Time%total
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%delta : ", Timer%Time%delta
            write(Test%outputUnit,"(*(g0))")   "Timer%period     : ", Timer%period
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%unit  : ", Timer%Time%unit
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_Timer_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getTimeSinceStart_1() result(assertion)

        use Constants_mod, only: IK, RK
        use System_mod, only: sleep
        implicit none
        real(RK), parameter :: seconds = 0.05_RK
        logical             :: assertion
        type(Timer_type)    :: Timer
        type(Err_type)      :: Err

        Timer = Timer_type(Err)
        assertion = .not. Err%occurred; if (.not. assertion) return
        call sleep(seconds=seconds,Err=Err)
        assertion = .not. Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Timer%total() > 0.9_RK * seconds

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%start: ", Timer%Count%start
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%rate : ", Timer%Count%rate
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%max  : ", Timer%Count%max
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%stop : ", Timer%Count%stop
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%total: ", Timer%Count%total
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%delta: ", Timer%Count%delta
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%start : ", Timer%Time%start
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%stop  : ", Timer%Time%stop
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%total : ", Timer%Time%total
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%delta : ", Timer%Time%delta
            write(Test%outputUnit,"(*(g0))")   "Timer%period     : ", Timer%period
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%unit  : ", Timer%Time%unit
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getTimeSinceStart_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getTimeSinceLastCall_1() result(assertion)

        use System_mod, only: sleep
        use Constants_mod, only: IK, RK
        implicit none
        logical             :: assertion
        real(RK), parameter :: seconds = 0.05_RK
        type(Timer_type)    :: Timer
        type(Err_type)      :: Err

        Timer = Timer_type(Err)
        assertion = .not. Err%occurred; if (.not. assertion) return
        call sleep(seconds=seconds,Err=Err)
        assertion = .not. Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Timer%delta() > 0.9_RK * seconds

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%start: ", Timer%Count%start
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%rate : ", Timer%Count%rate
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%max  : ", Timer%Count%max
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%stop : ", Timer%Count%stop
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%total: ", Timer%Count%total
            write(Test%outputUnit,"(*(g0))")   "Timer%Count%delta: ", Timer%Count%delta
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%start : ", Timer%Time%start
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%stop  : ", Timer%Time%stop
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%total : ", Timer%Time%total
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%delta : ", Timer%Time%delta
            write(Test%outputUnit,"(*(g0))")   "Timer%period     : ", Timer%period
            write(Test%outputUnit,"(*(g0))")   "Timer%Time%unit  : ", Timer%Time%unit
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getTimeSinceLastCall_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Timer_mod ! LCOV_EXCL_LINE