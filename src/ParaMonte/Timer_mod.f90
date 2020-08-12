!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module Timer_mod

    use, intrinsic :: iso_fortran_env, only: int64
    use Constants_mod, only: RK
    implicit none

    private
    public :: MODULE_NAME, Timer_type

    character(*), parameter :: MODULE_NAME = "@Timer_mod"

    type, private           :: Count_type
        integer(int64)      :: start, stop, total, delta, max
        real(RK)            :: rate
    end type Count_type

    type, private           :: Time_type
        real(RK)            :: start, stop, total, delta
        character(7)        :: unit = "seconds"
    end type Time_type

    type                    :: Timer_type
        type(Count_type)    :: Count
        type(Time_type)     :: Time
        real(RK)            :: period
    contains
        procedure, pass     :: tic => setTic
        procedure, pass     :: toc => setToc
        procedure, pass     :: delta => getTimeSinceLastCall
        procedure, pass     :: total => getTimeSinceStart
    end type Timer_type

    interface Timer_type
        module procedure    :: constructTimer
    end interface Timer_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructTimer(Err) result(Timer)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructTimer
#endif
        use Constants_mod, only: RK
        use Err_mod, only: Err_type
        implicit none
        type(Err_type), intent(out) :: Err
        type(Timer_type)            :: Timer
        character(*), parameter     :: PROCEDURE_NAME = "@constructTimer()"
        Err%occurred = .false.
        Err%msg = ""
        call system_clock( count=Timer%Count%start, count_rate=Timer%Count%rate, count_max=Timer%Count%max )
        if ( Timer%Count%start==-huge(0) .or. Timer%Count%rate==0._RK .or. Timer%Count%max==0 ) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred. There is no processor clock."
            return
        end if
        call Timer%tic()
    end function constructTimer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTic(Timer)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTic
#endif
        use Constants_mod, only: RK
        implicit none
        class(Timer_type), intent(inout)    :: Timer
        call system_clock( count=Timer%Count%start, count_rate=Timer%Count%rate, count_max=Timer%Count%max )
        Timer%period        = 1._RK / Timer%Count%rate
        Timer%Count%stop    = Timer%Count%start
        Timer%Count%total   = 0
        Timer%Count%delta   = 0
        Timer%Time%start    = real(Timer%Count%start,kind=RK) * Timer%period
        Timer%Time%stop     = Timer%Time%start
        Timer%Time%total    = 0._RK
        Timer%Time%delta    = 0._RK
    end subroutine setTic

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setToc(Timer)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setToc
#endif
        implicit none
        class(Timer_type), intent(inout)    :: Timer
        integer(int64)                      :: stopCount
        real(RK)                            :: stopTime
        call system_clock( count=stopCount )
        Timer%Count%delta   = stopCount - Timer%Count%stop
        Timer%Count%total   = stopCount - Timer%Count%start
        Timer%Count%stop    = stopCount
        stopTime            = real(stopCount,kind=RK) * Timer%period
        Timer%Time%delta    = stopTime - Timer%Time%stop
        Timer%Time%total    = stopTime - Timer%Time%start
        Timer%Time%stop     = stopTime
    end subroutine setToc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getTimeSinceLastCall(Timer) result(timeSinceLastCall)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeSinceLastCall
#endif
        implicit none
        class(Timer_type), intent(inout)    :: Timer
        real(RK)                            :: timeSinceLastCall
        call Timer%toc()
        timeSinceLastCall = Timer%Time%delta
    end function getTimeSinceLastCall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getTimeSinceStart(Timer) result(timeSinceStart)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeSinceStart
#endif
        implicit none
        class(Timer_type), intent(inout)    :: Timer
        real(RK)                            :: timeSinceStart
        call Timer%toc()
        timeSinceStart = Timer%Time%delta
    end function getTimeSinceStart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Timer_mod