!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module TimerCPU_mod

    use Err_mod, only: Err_type
    use Constants_mod, only: RK
    implicit none

    private
    public :: MODULE_NAME, TimerCPU_type

    character(*), parameter :: MODULE_NAME = "@TimerCPU_mod"

    type, private           :: Time_type
        real(RK)            :: start, stop, delta
        character(7)        :: unit = "seconds"
    end type Time_type

    type                    :: TimerCPU_type
        type(Time_type)     :: Time
        type(Err_type)      :: Err
    contains
        procedure, pass     :: tic => setTicCPU
        procedure, pass     :: toc => setTocCPU
    end type TimerCPU_type

    interface TimerCPU_type
        module procedure    :: constructTimerCPU
    end interface TimerCPU_type

!******************************************************************************************************************************************
!******************************************************************************************************************************************

contains

!******************************************************************************************************************************************
!******************************************************************************************************************************************

    function constructTimerCPU() result(TimerCPU)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructTimerCPU
#endif
        use Constants_mod, only: RK
        implicit none
        type(TimerCPU_type) :: TimerCPU
        character(*), parameter :: PROCEDURE_NAME = "@constructTimerCPU()"
        TimerCPU%Err%occurred = .false.
        TimerCPU%Err%msg = ""
        call cpu_time( time=TimerCPU%Time%start )
        if ( TimerCPU%Time%start<0 ) then
            TimerCPU%Err%occurred = .true.
            TimerCPU%Err%msg = PROCEDURE_NAME // ": There is no processor clock."
            return
        end if
        call TimerCPU%tic()
    end function constructTimerCPU

!******************************************************************************************************************************************
!******************************************************************************************************************************************

  subroutine setTicCPU(TimerCPU)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTicCPU
#endif
    implicit none
    class(TimerCPU_type), intent(inout)    :: TimerCPU
    call cpu_time( time=TimerCPU%Time%start )
  end subroutine setTicCPU

!******************************************************************************************************************************************
!******************************************************************************************************************************************

  subroutine setTocCPU(TimerCPU)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTocCPU
#endif
    implicit none
    class(TimerCPU_type), intent(inout)    :: TimerCPU
    call cpu_time( time=TimerCPU%Time%stop )
    TimerCPU%Time%delta = TimerCPU%Time%stop - TimerCPU%Time%start
  end subroutine setTocCPU

!******************************************************************************************************************************************
!******************************************************************************************************************************************

end module TimerCPU_mod