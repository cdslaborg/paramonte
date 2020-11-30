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

!>  \brief This module contains tests of the module [FileList_mod](@ref filelist_mod).
!>  @author Amir Shahmoradi

module Test_FileList_mod

    use FileList_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_FileList

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_FileList()

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getFileList_1, "test_getFileList_1")
        call Test%run(test_constructFileList_1, "test_constructFileList_1")
        call Test%run(test_constructFileList_2, "test_constructFileList_2")
        call Test%run(test_constructFileList_3, "test_constructFileList_3")
        call Test%run(test_constructFileList_4, "test_constructFileList_4")
        call Test%run(test_constructFileList_5, "test_constructFileList_5")
        call Test%run(test_constructFileList_6, "test_constructFileList_6")
        call Test%run(test_constructFileList_7, "test_constructFileList_7")
        call Test%run(test_constructFileList_8, "test_constructFileList_8")
        call Test%finalize()

    end subroutine test_FileList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileList_1() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type()
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. FileList%searchStr == ""
        assertion = assertion .and. FileList%orderStr == ""
        assertion = assertion .and. FileList%excludeStr == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileList_2() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type(orderStr="name")
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. FileList%searchStr == ""
        assertion = assertion .and. FileList%orderStr == "name"
        assertion = assertion .and. FileList%excludeStr == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileList_3() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type(orderStr="date")
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. FileList%searchStr == ""
        assertion = assertion .and. FileList%orderStr == "date"
        assertion = assertion .and. FileList%excludeStr == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileList_4() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type(excludeStr = "ParaMonte")
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. FileList%searchStr == ""
        assertion = assertion .and. FileList%orderStr == ""
        assertion = assertion .and. FileList%excludeStr == "ParaMonte"

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFileList_5() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type(searchStr = "ParaMonte")
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        assertion = assertion .and. FileList%searchStr == "ParaMonte"
        assertion = assertion .and. FileList%orderStr == ""
        assertion = assertion .and. FileList%excludeStr == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input argument `search` is optional and can be dropped, in which case, `searchStr` attribute becomes empty string `""`.
    function test_constructFileList_6() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type()
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        assertion = assertion .and. FileList%searchStr == ""
        assertion = assertion .and. FileList%orderStr == ""
        assertion = assertion .and. FileList%excludeStr == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The listing `order` must be either `"name"` or `"date"`. Anything else must yield error.
    function test_constructFileList_7() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        FileList = FileList_type(orderStr = "datename")
        assertion = FileList%Err%occurred
        if (.not. assertion) then
            ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
            ! LCOV_EXCL_STOP
        end if

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The constructor must be able to receive an optional `OS` argument.
    function test_constructFileList_8() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        use System_mod, only: OS_type
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i
        type(OS_type)       :: OS

        assertion = .true.

        call OS%query()
        FileList = FileList_type(OS = OS)
        assertion = .not. FileList%Err%occurred
        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "searchStr  : '", FileList%searchStr, "'"
            write(Test%outputUnit,"(*(g0))")   "orderStr   : '", FileList%orderStr, "'"
            write(Test%outputUnit,"(*(g0))")   "excludeStr : '", FileList%excludeStr, "'"
            do i = 1,FileList%count
                write(Test%outputUnit,"(*(g0))")   "FileList%File(" // num2str(i) // ") : ", FileList%File(i)%record
            end do
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFileList_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The [getFileList](@ref FileList_mod::getfilelist) procedure must be able to automatically set the missing optional arguments.
    function test_getFileList_1() result(assertion)
        use Constants_mod, only: IK
        use String_mod, only: num2str
        use System_mod, only: OS_type
        implicit none
        logical :: assertion
        type(FileList_type) :: FileList
        integer(IK)         :: i

        assertion = .true.

        call FileList%get(FileList%File, FileList%Err)
        assertion = .not. FileList%Err%occurred

        if (.not. assertion) then
        ! LCOV_EXCL_START
            if (Test%isDebugMode) then
                write(Test%outputUnit,"(*(g0))") "FileList%Err%occurred = ", FileList%Err%occurred
            end if
            return
        end if
        ! LCOV_EXCL_STOP

    end function test_getFileList_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_FileList_mod ! LCOV_EXCL_LINE