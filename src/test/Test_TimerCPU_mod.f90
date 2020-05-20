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

module Test_TimerCPU_mod

    use TimerCPU_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_TimerCPU

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_TimerCPU()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_TimerCPU_type()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif

        
    end subroutine test_TimerCPU
    
!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_TimerCPU_type()

        use System_mod, only: sleep
        use Constants_mod, only: RK
        implicit none
        type(TimerCPU_type) :: TimerCPU


        call Test%testing("TimerCPU_type")

        TimerCPU = TimerCPU_type()
        call Test%checkForErr(TimerCPU%Err)
        call sleep(seconds=0.1_RK,Err=TimerCPU%Err)
        call Test%checkForErr(TimerCPU%Err)
        call TimerCPU%toc()

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "TimerCPU%Time%start : ", TimerCPU%Time%start
            write(Test%outputUnit,"(*(g0))")   "TimerCPU%Time%stop  : ", TimerCPU%Time%stop
            write(Test%outputUnit,"(*(g0))")   "TimerCPU%Time%delta : ", TimerCPU%Time%delta
            write(Test%outputUnit,"(*(g0))")   "TimerCPU%Time%unit  : ", TimerCPU%Time%unit
            write(Test%outputUnit,"(*(g0))")
        end if
        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_TimerCPU_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_TimerCPU_mod