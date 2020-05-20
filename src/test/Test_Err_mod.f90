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

module Test_Err_mod

    use Err_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Err

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Err()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        if (Test%Image%isFirst) call test_note()
        if (Test%Image%isFirst) call test_warn()
        !if (Test%Image%isFirst) call test_abort()
        call Test%finalize()

    end subroutine test_Err

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_note()

        use Decoration_mod, only: TAB
        implicit none
        character(:), allocatable   :: prefix, msg
        prefix  = TAB//TAB//"ParaMonte"
        msg     = "What does a fish know about the water in which it swims all its life?    Albert Einstein\n" // &
                  "Everything should be made as simple as possible, but not simpler.    Albert Einstein\n" // &
                  "The absence of evidence is not evidence for absence.    Carl Sagan\n" // &
                  "If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton\n" // &
                  "I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        call Test%testing("@note()")

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)") "prefix: ", prefix
            write(Test%outputUnit,"(2A)") "msg   : ", msg
            write(Test%outputUnit,"(2A)")
            call note(msg,prefix,"\n")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_note

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_warn()

        use Decoration_mod, only: TAB
        implicit none
        character(:), allocatable   :: prefix, msg
        prefix  = TAB//TAB//"ParaMonte"
        msg     = "What does a fish know about the water in which it swims all its life?    Albert Einstein\n" // &
                  "Everything should be made as simple as possible, but not simpler.    Albert Einstein\n" // &
                  "The absence of evidence is not evidence for absence.    Carl Sagan\n" // &
                  "If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton\n" // &
                  "I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        call Test%testing("@warn()")

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)") "prefix: ", prefix
            write(Test%outputUnit,"(2A)") "msg   : ", msg
            write(Test%outputUnit,"(2A)")
            call warn(msg,prefix,"\n")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_warn

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_abort()

        use Decoration_mod, only: TAB
        implicit none
        character(:), allocatable   :: prefix, msg
        prefix  = TAB//TAB//"ParaMonte"
        msg     = "What does a fish know about the water in which it swims all its life?    Albert Einstein\n" // &
                  "Everything should be made as simple as possible, but not simpler.    Albert Einstein\n" // &
                  "The absence of evidence is not evidence for absence.    Carl Sagan\n" // &
                  "If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton\n" // &
                  "I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        call Test%testing("@abort()")

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(2A)")
            write(Test%outputUnit,"(2A)") "prefix: ", prefix
            write(Test%outputUnit,"(2A)") "msg   : ", msg
            write(Test%outputUnit,"(2A)")
            !call abort(msg,prefix,"\n")
        end if

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_abort

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Err_mod
