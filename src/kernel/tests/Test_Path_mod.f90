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

!>  \brief This module contains tests of the module [Path_mod](@ref path_mod).
!>  @author Amir Shahmoradi

module Test_Path_mod

    use Path_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Path

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Path()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_winify, "test_winify")
        call Test%run(test_linify_1, "test_linify_1")
        call Test%run(test_linify_2, "test_linify_2")
        call Test%run(test_constructPath, "test_constructPath")
        call Test%run(test_getDirNameExt_1, "test_getDirNameExt_1")
        call Test%run(test_getDirNameExt_2, "test_getDirNameExt_2")
        call Test%run(test_getDirFullName_1, "test_getDirFullName_1")
        call Test%run(test_getDirFullName_2, "test_getDirFullName_2")
        call Test%finalize()

    end subroutine test_Path

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPath() result(assertion)

        !use String_mod, only: num2str
        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        assertion = Path%original == "./temp\ Folder/\{inside\}\/"
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original  :", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%modified  :", Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%slashOS   :", Path%slashOS
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir       :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name      :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext       :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_constructPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_winify() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = "./temp\ Folder/\{inside\}/"
        call Path%winify(Path%original,Path%modified,Path%Err)
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        assertion = Path%modified == '".\temp Folder\{inside}\"'
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") Path%original
            write(Test%outputUnit,"(*(g0,:,' '))") Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_winify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_linify_1() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = '".\temp Folder\{inside}\"'
        call Path%linify(Path%original,Path%modified)

        assertion = Path%modified == "./temp\ Folder/\{inside\}/"

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original :", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%modified :", Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%slashOS  :", Path%slashOS
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir      :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name     :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext      :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_linify_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_linify_2() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = '".\temp Folder\{inside}\-"'
        call Path%linify(Path%original,Path%modified)

        assertion = Path%modified == "./temp\ Folder/\{inside\}/-"

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original:", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%modified:", Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%slashOS :", Path%slashOS
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir       :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name      :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext       :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if


    end function test_linify_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDirNameExt_1() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        assertion = .true.

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = ".\temp Folder\{inside}\"
        call Path%getDirNameExt(Path%original,"\",Path%dir,Path%name,Path%ext)

        assertion = assertion .and. Path%dir == ".\temp Folder\{inside}\"
        assertion = assertion .and. Path%name == ""
        assertion = assertion .and. Path%ext == ""

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getDirNameExt_1

    function test_getDirNameExt_2() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        assertion = .true.

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = ".\temp Folder\{inside}\Temp.tXt"
        call Path%getDirNameExt(Path%original,"\",Path%dir,Path%name,Path%ext)

        assertion = assertion .and. Path%dir == ".\temp Folder\{inside}\"
        assertion = assertion .and. Path%name == "Temp"
        assertion = assertion .and. Path%ext == ".tXt"

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getDirNameExt_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDirFullName_1() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        assertion = .true.

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = ".\temp Folder\{inside}\-"
        call Path%getDirFullName(Path%original,"\",Path%dir,Path%name)

        assertion = assertion .and. Path%dir == ".\temp Folder\{inside}\"
        assertion = assertion .and. Path%name == "-"

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getDirFullName_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDirFullName_2() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        assertion = .true.

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = "-"
        call Path%getNameExt(Path%original,Path%name,Path%ext)

        assertion = assertion .and. Path%name == "-"
        assertion = assertion .and. Path%ext == ""

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getDirFullName_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Path_mod
