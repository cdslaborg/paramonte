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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains classes and procedures relevant to wall-time timing.
!>  \author Amir Shahmoradi

module Timer_mod

    use, intrinsic :: iso_fortran_env, only: int64
    use Constants_mod, only: RK
    implicit none

    ! Doxygen does not comprehend the following statements. Therefore, these are commented out.
    !private
    !public :: MODULE_NAME, Timer_type

    character(*), parameter :: MODULE_NAME = "@Timer_mod"

    type                    :: Count_type
        integer(int64)      :: start            !< The first processor clock count.
        integer(int64)      :: stop             !< The last fetched processor clock count.
        integer(int64)      :: total            !< The total processor clock count since start.
        integer(int64)      :: delta            !< The total processor clock count since the last measurement.
        integer(int64)      :: max              !< The maximum value that the processor count may take, or zero if there is no clock.
        real(RK)            :: rate             !< The number of clock counts per second, or zero if there is no clock.
    end type Count_type

    type                    :: Time_type
        real(RK)            :: start            !< The start time in seconds.
        real(RK)            :: stop             !< The stop time in seconds.
        real(RK)            :: total            !< The total time in seconds since the start.
        real(RK)            :: delta            !< The total time in seconds since the last timing.
        character(7)        :: unit = "seconds" !< The unit of time.
    end type Time_type

    !> The `Timer_type` class, containing method for setting up a wall-time timer.
    type                    :: Timer_type
        type(Count_type)    :: Count            !< An object of type [Count_type](@ref count_type) containing information about the processor clock count.
        type(Time_type)     :: Time             !< An object of type [Time_type](@ref time_type) containing information about the processor time.
        real(RK)            :: period           !< The time between the processor clock tics in seconds.
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

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is the constructor of the class [Timer_type](@ref timer_type).
    !> Before returning the object, this function also calls the [tic()](@ref settic)
    !> method of the `Timer_type` object to reset the timer.
    !>
    !> \param[out]  Err :   An object of class [Err_type](@ref err_mod::err_type) indicating
    !>                      the occurrence of error during the object construction.
    !>
    !> \return
    !> `Timer` : An object of class [Timer_type](@ref timer_type).
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    function constructTimer(Err) result(Timer)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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
            ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred. There is no processor clock."
            return
            ! LCOV_EXCL_STOP
        end if
        call Timer%tic()
    end function constructTimer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the class [Timer_type](@ref timer_type).
    !> Reset the timer object and return.
    !>
    !> \param[inout]    Timer   :   An object of class [Timer_type](@ref timer_type) whose components are manipulated by this method.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    subroutine setTic(Timer)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> This procedure is a method of the class [Timer_type](@ref timer_type).
    !> Mark the timer and compute the time spent since last [toc()](@ref settoc) method call and return.
    !>
    !> \param[inout]    Timer   :   An object of class [Timer_type](@ref timer_type) whose components are manipulated by this method.
    !>
    !> \remark
    !> Specifically, this method will set/update the following components of the object of type `Timer_type`:
    !> + `Timer%Count%delta`    : The total counts since the last `toc()` call.
    !> + `Timer%Count%total`    : The total counts since the object creation or since the last `tic()` call.
    !> + `Timer%Count%stop`     : The current count as inferred from the Fortran intrinsic procedure `system_clock()`.
    !> + `Timer%Timer%delta`    : The total time in seconds since the last `toc()` call.
    !> + `Timer%Timer%total`    : The total time in seconds since the object creation or since the last `tic()` call.
    !> + `Timer%Timer%stop`     : The current time in seconds as inferred from the Fortran intrinsic procedure `system_clock()`.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    subroutine setToc(Timer)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> This procedure is a method of the class [Timer_type](@ref timer_type).
    !> Report the time spent in seconds since the last timing.
    !>
    !> \param[inout]    Timer   :   An object of class [Timer_type](@ref timer_type) whose components are manipulated by this method.
    !>
    !> \return
    !> `timeSinceLastCall`      :   The time spent in seconds since the last timing.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    function getTimeSinceLastCall(Timer) result(timeSinceLastCall)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeSinceLastCall
#endif
        implicit none
        class(Timer_type), intent(inout)    :: Timer
        real(RK)                            :: timeSinceLastCall
        call Timer%toc()
        timeSinceLastCall = Timer%Time%delta
    end function getTimeSinceLastCall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the class [Timer_type](@ref timer_type).
    !> Report the time spent in seconds since the start of the timing.
    !>
    !> \param[inout]    Timer   :   An object of class [Timer_type](@ref timer_type) whose components are manipulated by this method.
    !>
    !> \return
    !> `timeSinceStart`         :   The time spent in seconds since the start of the timing.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    function getTimeSinceStart(Timer) result(timeSinceStart)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getTimeSinceStart
#endif
        implicit none
        class(Timer_type), intent(inout)    :: Timer
        real(RK)                            :: timeSinceStart
        call Timer%toc()
        timeSinceStart = Timer%Time%delta
    end function getTimeSinceStart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Timer_mod ! LCOV_EXCL_LINE