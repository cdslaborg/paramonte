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

!>  \brief This module contains tests of the module [DateTime_mod](@ref datetime_mod).
!>  @author Amir Shahmoradi

module Test_DateTime_mod

    use DateTime_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_DateTime

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_DateTime()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_queryDateTime, "test_queryDateTime")
        call Test%run(test_getNiceDateTime, "test_getNiceDateTime")
        call Test%finalize()

    end subroutine test_DateTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_queryDateTime() result(assertion)

        implicit none
        logical             :: assertion
        type(DateTime_type) :: DateTime

        assertion = .true.

        call DateTime%query()

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "           date:", DateTime%date
            write(Test%outputUnit,"(*(g0,:,' '))") "           time:", DateTime%time
            write(Test%outputUnit,"(*(g0,:,' '))") "           zone:", DateTime%zone
            write(Test%outputUnit,"(*(g0,:,' '))") "        century:", DateTime%century
            write(Test%outputUnit,"(*(g0,:,' '))") "           year:", DateTime%year
            write(Test%outputUnit,"(*(g0,:,' '))") "          month:", DateTime%month
            write(Test%outputUnit,"(*(g0,:,' '))") "            day:", DateTime%day
            write(Test%outputUnit,"(*(g0,:,' '))") "           hour:", DateTime%hour
            write(Test%outputUnit,"(*(g0,:,' '))") "         minute:", DateTime%minute
            write(Test%outputUnit,"(*(g0,:,' '))") "         second:", DateTime%second
            write(Test%outputUnit,"(*(g0,:,' '))") "    millisecond:", DateTime%millisecond
            write(Test%outputUnit,"(*(g0,:,' '))") "         values:", DateTime%Values
            write(Test%outputUnit,"(*(g0,:,' '))") "fancyStyleBasic:", DateTime%fancyStyleBasic
            write(Test%outputUnit,"(*(g0,:,' '))") "     fancyStyle:", DateTime%fancyStyle
            write(Test%outputUnit,"(*(g0,:,' '))") 
        end if

    end function test_queryDateTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Test the equivalence of [getNiceDateTime()](@ref datetime_mod::getnicedatetime) with
    ! the output of the `fancyStyleBasic()` method of [DateTime_type](@ref datetime_mod::datetime_type).
    function test_getNiceDateTime() result(assertion)

        implicit none
        logical                     :: assertion
        type(DateTime_type)         :: DateTime
        character(:), allocatable   :: niceDateTime

        call DateTime%query()
        niceDateTime = getNiceDateTime()
        assertion = niceDateTime == DateTime%fancyStyleBasic

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "       getNiceDateTime():", niceDateTime
            write(Test%outputUnit,"(*(g0,:,' '))") "DateTime%fancyStyleBasic:", DateTime%fancyStyleBasic
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

        assertion = .true.

    end function test_getNiceDateTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_DateTime_mod
