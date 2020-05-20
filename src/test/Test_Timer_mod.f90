!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module Test_Timer_mod

    use Timer_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Timer

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Timer()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_Timer_type()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif


    end subroutine test_Timer

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Timer_type()
        
        use System_mod, only: sleep
        use Constants_mod, only: RK
        implicit none
        type(Timer_type)    :: Timer
        type(Err_type)      :: Err


        if (Test%Image%isFirst) call Test%testing("Timer_type")

        Timer = Timer_type(Err)
        call Test%checkForErr(Err)
        call sleep(seconds=0.1_RK,Err=Err)
        call Test%checkForErr(Err)
        call Timer%toc()

        if (Test%isDebugMode .and. Test%Image%isFirst) then
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
        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_Timer_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Timer_mod