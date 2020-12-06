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

!>  \brief This module contains tests of the module [Err_mod](@ref err_mod).
!>  \author Amir Shahmoradi

module Test_Err_mod

    use Err_mod
    use Decoration_mod, only: TAB
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Err

    type(Test_type) :: Test

    character(*), parameter :: mc_prefix =  TAB//TAB//"ParaMonte"
    character(*), parameter :: mc_msg    =  "What does a fish know about the water in which it swims all its life?    Albert Einstein\n" // &
                                            "Everything should be made as simple as possible, but not simpler.    Albert Einstein\n" // &
                                            "The absence of evidence is not evidence for absence.    Carl Sagan\n" // &
                                            "If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton\n" // &
                                            "I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Err()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_note_1, "test_note_1")
        call Test%run(test_note_2, "test_note_2")
        call Test%run(test_warn_1, "test_warn_1")
        call Test%run(test_warn_2, "test_warn_2")
        call Test%run(test_abort_1, "test_abort_1")
        call Test%run(test_abort_2, "test_abort_2")
        call Test%finalize()

    end subroutine test_Err

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_note_1() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        use String_mod, only: num2str

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 5_IK
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = "        ParaMonte - NOTE: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        OutputList_ref(2)%record = "        ParaMonte - NOTE: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        OutputList_ref(3)%record = "        ParaMonte - NOTE: The absence of evidence is not evidence for absence.    Carl Sagan"
        OutputList_ref(4)%record = "        ParaMonte - NOTE: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        OutputList_ref(5)%record = "        ParaMonte - NOTE: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"/Test_Err_mod@test_note_1."//num2str(Test%Image%id)//".out", status = "replace")

        call note   ( msg = mc_msg &
                    , prefix  = mc_prefix &
                    , newline  = "\n" &
                    , outputUnit = fileUnit &
                    , marginTop = 0_IK &
                    , marginBot = 0_IK &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)

            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            assertion = iostat == 0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            OutputList(i)%record = trim(OutputList(i)%record)

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_note_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_note_2() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK, NLC
        use String_mod, only: num2str, replaceStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 7_IK
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = ""
        OutputList_ref(2)%record = " - NOTE: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        OutputList_ref(3)%record = " - NOTE: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        OutputList_ref(4)%record = " - NOTE: The absence of evidence is not evidence for absence.    Carl Sagan"
        OutputList_ref(5)%record = " - NOTE: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        OutputList_ref(6)%record = " - NOTE: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"
        OutputList_ref(7)%record = ""

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"/Test_Err_mod@test_note_2."//num2str(Test%Image%id)//".out", status = "replace")

        call note   ( msg = replaceStr(mc_msg, "\n", NLC) &
                    , outputUnit = fileUnit &
                    , newline  = NLC &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)

            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            assertion = iostat == 0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            OutputList(i)%record = trim(OutputList(i)%record)

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_note_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_warn_1() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        use String_mod, only: num2str

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 5_IK
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = "        ParaMonte - WARNING: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        OutputList_ref(2)%record = "        ParaMonte - WARNING: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        OutputList_ref(3)%record = "        ParaMonte - WARNING: The absence of evidence is not evidence for absence.    Carl Sagan"
        OutputList_ref(4)%record = "        ParaMonte - WARNING: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        OutputList_ref(5)%record = "        ParaMonte - WARNING: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"/Test_Err_mod@test_warn_1."//num2str(Test%Image%id)//".out", status = "replace")

        call warn   ( msg = mc_msg &
                    , prefix  = mc_prefix &
                    , newline  = "\n" &
                    , outputUnit = fileUnit &
                    , marginTop = 0_IK &
                    , marginBot = 0_IK &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)

            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            assertion = iostat == 0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            OutputList(i)%record = trim(OutputList(i)%record)

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_warn_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_warn_2() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK, NLC
        use String_mod, only: num2str, replaceStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 7_IK
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)

        assertion = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = ""
        OutputList_ref(2)%record = " - WARNING: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        OutputList_ref(3)%record = " - WARNING: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        OutputList_ref(4)%record = " - WARNING: The absence of evidence is not evidence for absence.    Carl Sagan"
        OutputList_ref(5)%record = " - WARNING: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        OutputList_ref(6)%record = " - WARNING: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"
        OutputList_ref(7)%record = ""

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"/Test_Err_mod@test_warn_2."//num2str(Test%Image%id)//".out", status = "replace")

        call warn   ( msg = replaceStr(mc_msg, "\n", NLC) &
                    , outputUnit = fileUnit &
                    , newline  = NLC &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)

            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            assertion = iostat == 0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            OutputList(i)%record = trim(OutputList(i)%record)

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_warn_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_abort_1() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        use String_mod, only: num2str

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 6_IK
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)
        type(Err_type) :: Err

        assertion = .true.
        mv_isTestingMode = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = ""
        OutputList_ref(2)%record = "        ParaMonte - FATAL: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        OutputList_ref(3)%record = "        ParaMonte - FATAL: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        OutputList_ref(4)%record = "        ParaMonte - FATAL: The absence of evidence is not evidence for absence.    Carl Sagan"
        OutputList_ref(5)%record = "        ParaMonte - FATAL: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        OutputList_ref(6)%record = "        ParaMonte - FATAL: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"/Test_Err_mod@test_abort_1."//num2str(Test%Image%id)//".out", status = "replace")

        Err%msg = mc_msg
        call abort  ( Err = Err &
                    , prefix  = mc_prefix &
                    , newline  = "\n" &
                    , outputUnit = fileUnit &
                    , returnEnabled = .true. &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)

            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            assertion = iostat == 0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            OutputList(i)%record = trim(OutputList(i)%record)

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = ", OutputList_ref(i)%record
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = ", OutputList(i)%record
                write(Test%outputUnit,"(*(g0))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_abort_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the effects of an input non-null error code `Err%stat`.
    !> Test the effects of missing arguments `prefix`, `returnEnabled`, and `newline`.
    function test_abort_2() result(assertion)

        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK, NLC
        use String_mod, only: num2str, replaceStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 7_IK
        type(CharVec_type), allocatable :: OutputList_ref(:)
        type(CharVec_type), allocatable :: OutputList(:)
        type(Err_type) :: Err

        assertion = .true.
        mv_isTestingMode = .true.

        if (allocated(OutputList)) deallocate(OutputList); allocate(OutputList(NLINE))
        if (allocated(OutputList_ref)) deallocate(OutputList_ref); allocate(OutputList_ref(NLINE))

        OutputList_ref(1)%record = ""
        OutputList_ref(2)%record = " - FATAL: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        OutputList_ref(3)%record = " - FATAL: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        OutputList_ref(4)%record = " - FATAL: The absence of evidence is not evidence for absence.    Carl Sagan"
        OutputList_ref(5)%record = " - FATAL: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        OutputList_ref(6)%record = " - FATAL: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"
        OutputList_ref(7)%record = " - FATAL: Error Code: 123."

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = Test%outDir//"/Test_Err_mod@test_abort_2."//num2str(Test%Image%id)//".out", status = "replace")

        Err%msg = replaceStr(mc_msg, "\n", NLC)
        Err%stat = 123_IK
        call abort  ( Err = Err &
                    , outputUnit = fileUnit &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(OutputList(i)%record)) deallocate(OutputList(i)%record)
            allocate(character(132) :: OutputList(i)%record)

            read(fileUnit,"(A132)", iostat = iostat) OutputList(i)%record
            assertion = iostat == 0_IK
            if (.not. assertion) return ! LCOV_EXCL_LINE

            OutputList(i)%record = trim(OutputList(i)%record)

            assertionCurrent = OutputList(i)%record == OutputList_ref(i)%record
            assertion = assertion .and. assertionCurrent

            if (Test%isDebugMode .and. .not. assertionCurrent) then
            ! LCOV_EXCL_START
                write(Test%outputUnit,"(*(g0))")
                write(Test%outputUnit,"(*(g0))") "OutputList_ref(",num2str(i),")%record = '", OutputList_ref(i)%record, "'"
                write(Test%outputUnit,"(*(g0))") "OutputList    (",num2str(i),")%record = '", OutputList(i)%record, "'"
                write(Test%outputUnit,"(*(g0))")
            end if
            ! LCOV_EXCL_STOP

        end do

    end function test_abort_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Err_mod ! LCOV_EXCL_LINE