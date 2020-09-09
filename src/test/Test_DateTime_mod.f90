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
