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

module Test_DateTime_mod

    use DateTime_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_DateTime

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_DateTime()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        if (Test%Image%isFirst) call test_DateTime_type()
        call Test%finalize()

    end subroutine test_DateTime

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_DateTime_type()

        implicit none
        type(DateTime_type) :: DateTime

        if (Test%Image%isFirst) call Test%testing("DateTime_type")
        call DateTime%query()

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)")           "Date:             ", DateTime%date
            write(Test%outputUnit,"(2A)")           "time:             ", DateTime%time
            write(Test%outputUnit,"(2A)")           "zone:             ", DateTime%zone
            write(Test%outputUnit,"(2A)")           "century:          ", DateTime%century
            write(Test%outputUnit,"(2A)")           "year:             ", DateTime%year
            write(Test%outputUnit,"(2A)")           "month:            ", DateTime%month
            write(Test%outputUnit,"(2A)")           "day:              ", DateTime%day
            write(Test%outputUnit,"(2A)")           "hour:             ", DateTime%hour
            write(Test%outputUnit,"(2A)")           "minute:           ", DateTime%minute
            write(Test%outputUnit,"(2A)")           "second:           ", DateTime%second
            write(Test%outputUnit,"(2A)")           "millisecond:      ", DateTime%millisecond
            write(Test%outputUnit,"(A,*(g0,' '))")  "values:           ", DateTime%Values
            write(Test%outputUnit,"(2A)")           "fancyStyleBasic:  ", DateTime%fancyStyleBasic
            write(Test%outputUnit,"(2A)")           "fancyStyle:       ", DateTime%fancyStyle
            write(Test%outputUnit,"(2A)")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

        call Test%testing("equivalence of getNiceDateTime() with DateTime_type%fancyStyleBasic()")
        call DateTime%query()
        Test%assertion = getNiceDateTime() == DateTime%fancyStyleBasic
        call Test%verify()

    end subroutine test_DateTime_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_DateTime_mod
