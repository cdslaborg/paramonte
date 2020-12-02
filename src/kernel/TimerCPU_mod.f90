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

!>  \brief This module contains classes and procedures relevant to CPU-time timing.
!>  @author Amir Shahmoradi

module TimerCPU_mod

    use Err_mod, only: Err_type
    use Constants_mod, only: RK
    implicit none

    ! Doxygen does not comprehend the following statements. Therefore, are commented out.
    !private
    !public :: MODULE_NAME, TimerCPU_type

    character(*), parameter :: MODULE_NAME = "@TimerCPU_mod"

    type                    :: Time_type
        real(RK)            :: start            !< The CPU start time.
        real(RK)            :: stop             !< The CPU stop time.
        real(RK)            :: delta            !< The CPU time spent since the last clock call in seconds.
        real(RK)            :: total            !< The total CPU time spent between start and stop in seconds.
        character(7)        :: unit = "seconds" !< The time unit.
    end type Time_type

    !> The `TimerCPU_type` class, containing method for setting up a CPU-time timer.
    type                    :: TimerCPU_type
        type(Time_type)     :: Time !< An object of type [Time_type](@ref time_type) containing information about the processor clock count.
        type(Err_type)      :: Err  !< An object of type [Err_type](@ref err_mod::err_type) containing error-handling information.
    contains
        procedure, pass     :: tic => setTicCPU
        procedure, pass     :: toc => setTocCPU
    end type TimerCPU_type

    interface TimerCPU_type
        module procedure    :: constructTimerCPU
    end interface TimerCPU_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is the constructor of the class [TimerCPU_type](@ref timercpu_type).
    !> Before returning the object, this function also calls the [tic()](@ref setticcpu)
    !> method of the object of class [TimerCPU_type](@ref timercpu_type) to reset the timer.
    !>
    !> \return
    !> `Timer` : An object of class [TimerCPU_type](@ref timercpu_type).
    !>
    !> \warning
    !> Upon return, be sure to check the value of `TimerCPU_type%%Err%%occurred` for the occurrence of any potential error.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    function constructTimerCPU() result(TimerCPU)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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
        ! LCOV_EXCL_START
            TimerCPU%Err%occurred = .true.
            TimerCPU%Err%msg = PROCEDURE_NAME // ": There is no processor clock."
            return
        end if
        ! LCOV_EXCL_STOP
        call TimerCPU%tic()
    end function constructTimerCPU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the class [TimerCPU_type](@ref timercpu_type).
    !> Reset the timer object and return.
    !>
    !> \param[inout] TimerCPU : An object of class [TimerCPU_type](@ref timercpu_type) whose components are manipulated by this method.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    subroutine setTicCPU(TimerCPU)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTicCPU
#endif
        implicit none
        class(TimerCPU_type), intent(inout) :: TimerCPU
        call cpu_time( time=TimerCPU%Time%start )
        TimerCPU%Time%delta = 0._RK
        TimerCPU%Time%total = 0._RK
        TimerCPU%Time%stop = TimerCPU%Time%start
    end subroutine setTicCPU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the class [TimerCPU_type](@ref timercpu_type).
    !> Mark the CPU timer and compute the CPU time spent since last [toc()](@ref settoccpu) method call and return.
    !>
    !> \param[inout] TimerCPU : An object of class [TimerCPU_type](@ref timercpu_type) whose components are manipulated by this method.
    !>
    !> \remark
    !> Specifically, this method will set/update the following components of the object of type `Timer_type`:
    !> + `TimerCPU%Timer%stop`     : The current CPU time in seconds as inferred from the Fortran intrinsic procedure `cpu_time()`.
    !> + `TimerCPU%Timer%delta`    : The total time in seconds since the last `toc()` call.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    subroutine setTocCPU(TimerCPU)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTocCPU
#endif
        use Constants_mod, only: RK
        implicit none
        class(TimerCPU_type), intent(inout) :: TimerCPU
        real(RK)                            :: dummy
        call cpu_time( time=TimerCPU%Time%stop )
        dummy = TimerCPU%Time%stop - TimerCPU%Time%start
        TimerCPU%Time%delta = dummy - TimerCPU%Time%total
        TimerCPU%Time%total = dummy
    end subroutine setTocCPU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module TimerCPU_mod ! LCOV_EXCL_LINE