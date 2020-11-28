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

!>  \brief This module contains classes and procedures for computing the date and time.
!>  @author Amir Shahmoradi

module DateTime_mod

    implicit none

    public
    private :: queryDateTime

    character(*), parameter :: MODULE_NAME = "@DateTime_mod"

    !> The datetime_type class
    type, public :: DateTime_type
        character(8)    :: date
        character(10)   :: time
        character(5)    :: zone
        integer         :: Values(8)
        character(2)    :: century
        character(4)    :: year
        character(2)    :: month
        character(2)    :: day
        character(2)    :: hour
        character(2)    :: minute
        character(2)    :: second
        character(3)    :: millisecond
        character(21)   :: fancyStyleBasic
        character(35)   :: fancyStyle
    contains
        !> Query date and time.
        procedure, pass :: query => queryDateTime
    end type DateTime_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Function to query the [DateTime_type](@ref datetime_type) instance to the current time.
    !>
    !> @param[in,out]   self : An object of [DateTime_type](@ref datetime_type) class.
    subroutine queryDateTime(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: queryDateTime
#endif
        implicit none
        class(DateTime_type), intent(inout) :: self
        call date_and_time( date    = self%date   &
                          , time    = self%time   &
                          , zone    = self%zone   &
                          , Values  = self%Values &
                          )
        self%century            = self%date(1:2)
        self%year               = self%date(1:4)
        self%month              = self%date(5:6)
        self%day                = self%date(7:8)
        self%hour               = self%time(1:2)
        self%minute             = self%time(3:4)
        self%second             = self%time(5:6)
        self%millisecond        = self%time(8:10)
        self%fancyStyleBasic    = self%year // "/" // &             ! 5 characters
                                  self%month // "/" // &            ! 3 characters
                                  self%day //" - " // &             ! 5 characters
                                  self%hour // ":" // &             ! 3 characters
                                  self%minute // ":" // &           ! 3 characters
                                  self%second                       ! 2 characters
        self%fancyStyle         = self%fancyStyleBasic // "." // &  ! basic+1 characters
                                  self%millisecond // " " // &      ! 4 characters
                                  self%zone // " UTC"               ! 9 characters
    end subroutine queryDateTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return date and time in a nice format.
    !>
    !> \return
    !> A character vector containing date and time in a nice format.
    !>
    !> \remark
    !> This is an impure function due to its dependence on `date_and_time()` Fortran intrinsic function.
    function getNiceDateTime() result(niceDateTime)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNiceDateTime
#endif
        implicit none
        character(len=21)                        :: niceDateTime
        character(10)                            :: time
        character(8)                             :: date
        call date_and_time(date,time)
        niceDateTime = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' - '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    end function getNiceDateTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module DateTime_mod ! LCOV_EXCL_LINE