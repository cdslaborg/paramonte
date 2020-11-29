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

!>  \brief This module contains tests of the module [FileContents_mod](@ref filecontents_mod).
!>  @author Amir Shahmoradi

module Test_FileContents_mod

    use FileContents_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_FileContents

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_FileContents()

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_constructFileContents_1, "test_constructFileContents_1")
        call Test%run(test_constructFileContents_2, "test_constructFileContents_2")
        call Test%finalize()

    end subroutine test_FileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileContents_1() result(assertion)
        use System_mod, only: RandomFileName_type, removeFile, OS_type
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RandomFileName
        type(FileContents_type)     :: FileContents
        type(OS_type)               :: OS
        integer                     :: fileUnit, i, iostat

        RandomFileName = RandomFileName_type(key = Test%outDir//"Test_FileContents_mod@test_constructFileContents")
        assertion = .not. RandomFileName%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "RandomFileName%Err%occurred = ", RandomFileName%Err%occurred
                write(Test%outputUnit,"(*(g0))") "RandomFileName%Err%stat = ", RandomFileName%Err%stat
                write(Test%outputUnit,"(*(g0))") "RandomFileName%Err%msg = ", RandomFileName%Err%msg
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        open(newunit = fileUnit, file = RandomFileName%path, status = "new", iostat = iostat)
        assertion = iostat == 0
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "iostat = ", iostat, " <= 0"
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        write(fileUnit,"(A)") "Testing FileContents_type..."
        write(fileUnit,"(A)") " "
        write(fileUnit,"(A)") ""
        write(fileUnit,"(A)")

        FileContents = FileContents_type(RandomFileName%path)
        assertion = .not. FileContents%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileContents%Err%occurred = ", FileContents%Err%occurred
                write(Test%outputUnit,"(*(g0))") "FileContents%Err%stat = ", FileContents%Err%stat
                write(Test%outputUnit,"(*(g0))") "FileContents%Err%msg = ", FileContents%Err%msg
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        call OS%query()
        assertion = .not. OS%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "OS%Err%occurred = ", OS%Err%occurred
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        call removeFile(RandomFileName%path,OS%Err)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "OS%Err%occurred = ", OS%Err%occurred
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        assertion = assertion .and. FileContents%numRecord==4
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileContents%numRecord = ", FileContents%numRecord, " /= ", 4
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        assertion = assertion .and. FileContents%Line(1)%record=="Testing FileContents_type..."
        assertion = assertion .and. FileContents%Line(2)%record==""
        assertion = assertion .and. FileContents%Line(3)%record==""
        assertion = assertion .and. FileContents%Line(4)%record==""

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "numRecord  : ", num2str(FileContents%numRecord)
            do i = 1,FileContents%numRecord
                write(Test%outputUnit,"(*(g0))") "FileContents%Line(" // num2str(i) // ")%record : '", FileContents%Line(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileContents_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileContents_2() result(assertion)

        use System_mod, only: RandomFileName_type, removeFile, OS_type
        use String_mod, only: num2str
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RandomFileName
        type(FileContents_type)     :: FileContents
        type(OS_type)               :: OS
        integer                     :: fileUnit, i, iostat
        logical                     :: exist

        RandomFileName = RandomFileName_type(key = Test%outDir//"Test_FileContents_mod@test_constructFileContents")
        assertion = .not. RandomFileName%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "RandomFileName%Err%occurred = ", RandomFileName%Err%occurred
                write(Test%outputUnit,"(*(g0))") "RandomFileName%Err%stat = ", RandomFileName%Err%stat
                write(Test%outputUnit,"(*(g0))") "RandomFileName%Err%msg = ", RandomFileName%Err%msg
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        open(newunit = fileUnit, file = RandomFileName%path, status = "new", iostat = iostat)
        assertion = iostat == 0
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "iostat = ", iostat, " <= 0"
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        write(fileUnit,"(A)") "Testing FileContents_type..."
        write(fileUnit,"(A)") " "
        write(fileUnit,"(A)") ""
        write(fileUnit,"(A)")

        FileContents = FileContents_type(RandomFileName%path, delEnabled = .true.)
        assertion = .not. FileContents%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileContents%Err%occurred = ", FileContents%Err%occurred
                write(Test%outputUnit,"(*(g0))") "FileContents%Err%stat = ", FileContents%Err%stat
                write(Test%outputUnit,"(*(g0))") "FileContents%Err%msg = ", FileContents%Err%msg
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        call OS%query()
        assertion = .not. OS%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "OS%Err%occurred = ", OS%Err%occurred
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        inquire(file = RandomFileName%path, exist = exist)
        assertion = .not. exist
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "delEnabled = ", .true.
                write(Test%outputUnit,"(*(g0))") "exist = ", exist
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        assertion = assertion .and. FileContents%numRecord==4
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileContents%numRecord = ", FileContents%numRecord, " /= ", 4
            end if
            ! LCOV_EXCL_STOP
            return
        end if

        assertion = assertion .and. FileContents%Line(1)%record=="Testing FileContents_type..."
        assertion = assertion .and. FileContents%Line(2)%record==""
        assertion = assertion .and. FileContents%Line(3)%record==""
        assertion = assertion .and. FileContents%Line(4)%record==""

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "numRecord  : ", num2str(FileContents%numRecord)
            do i = 1,FileContents%numRecord
                write(Test%outputUnit,"(*(g0))") "FileContents%Line(" // num2str(i) // ")%record : '", FileContents%Line(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileContents_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_FileContents_mod ! LCOV_EXCL_LINE