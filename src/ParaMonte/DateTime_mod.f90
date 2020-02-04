!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
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

module DateTime_mod

    implicit none
    
    public
    private :: query

    character(*), parameter :: MODULE_NAME = "@DateTime_mod"

    type :: DateTime_type
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
        procedure, pass :: query
    end type DateTime_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine query(DateTime)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: query
#endif
        implicit none
        class(DateTime_type), intent(inout)    :: DateTime
        call date_and_time( date    = DateTime%date   &
                          , time    = DateTime%time   &
                          , zone    = DateTime%zone   &
                          , Values  = DateTime%Values &
                          )
        DateTime%century            = DateTime%date(1:2)
        DateTime%year               = DateTime%date(1:4)
        DateTime%month              = DateTime%date(5:6)
        DateTime%day                = DateTime%date(7:8)
        DateTime%hour               = DateTime%time(1:2)
        DateTime%minute             = DateTime%time(3:4)
        DateTime%second             = DateTime%time(5:6)
        DateTime%millisecond        = DateTime%time(8:10)
        DateTime%fancyStyleBasic    = DateTime%year // "/" // &    ! 5 characters
                                      DateTime%month // "/" // &   ! 3 characters
                                      DateTime%day //" - " // &    ! 5 characters
                                      DateTime%hour // ":" // &    ! 3 characters
                                      DateTime%minute // ":" // &  ! 3 characters
                                      DateTime%second              ! 2 characters
        DateTime%fancyStyle         = DateTime%fancyStyleBasic // "." // & ! basic+1 characters
                                      DateTime%millisecond // " " // &     ! 4 characters
                                      DateTime%zone // " UTC"              ! 9 characters
    end subroutine query

!***********************************************************************************************************************************
!***********************************************************************************************************************************
   
    function getNiceDateTime()
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNiceDateTime
#endif
        implicit none
        character(len=21)                        :: getNiceDateTime
        character(10)                            :: time
        character(8)                             :: date
        call date_and_time(date,time)
        getNiceDateTime = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' - '//time(1:2)//':'//time(3:4)//':'//time(5:6)
    end function getNiceDateTime

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module DateTime_mod
