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
!>  \author Amir Shahmoradi

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
        call Test%run(test_isdir_1, "test_isdir_1")
        call Test%run(test_query_1, "test_query_1")
        call Test%run(test_query_2, "test_query_2")
        call Test%run(test_query_3, "test_query_3")
        call Test%run(test_mkdir_1, "test_mkdir_1")
        call Test%run(test_mkdir_2, "test_mkdir_2")
        call Test%run(test_mkdir_3, "test_mkdir_3")
        call Test%run(test_modify_1, "test_modify_1")
        call Test%run(test_winify_1, "test_winify_1")
        call Test%run(test_winify_2, "test_winify_2")
        call Test%run(test_winify_3, "test_winify_3")
        call Test%run(test_linify_1, "test_linify_1")
        call Test%run(test_linify_2, "test_linify_2")
        call Test%run(test_getNameExt_1, "test_getNameExt_1")
        call Test%run(test_getNameExt_2, "test_getNameExt_2")
        call Test%run(test_constructPath, "test_constructPath")
        call Test%run(test_getDirNameExt_1, "test_getDirNameExt_1")
        call Test%run(test_getDirNameExt_2, "test_getDirNameExt_2")
        call Test%run(test_getDirFullName_1, "test_getDirFullName_1")
        call Test%run(test_getDirFullName_2, "test_getDirFullName_2")
        call Test%finalize()

    end subroutine test_Path

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test if `isdir()` can successfully detect an existing directory.
    function test_isdir_1() result(assertion)
        use Constants_mod, only: RK
        implicit none
        logical         :: assertion
        assertion = isdir("../")
    end function test_isdir_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When the original path is not allocated, `guery()` must return an error message.
    function test_query_1() result(assertion)

        !use String_mod, only: num2str
        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        deallocate(Path%original)
        call Path%query()
        assertion = Path%Err%occurred

    end function test_query_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When the original path is allocated but is empty, `guery()` must return an error message.
    function test_query_2() result(assertion)

        !use String_mod, only: num2str
        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = ""
        call Path%query()
        assertion = Path%Err%occurred

    end function test_query_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When the optional OS is provided, the results must be the same as when it is not provided.
    function test_query_3() result(assertion)

        use System_mod, only: OS_type
        implicit none
        logical         :: assertion
        type(Path_type) :: Path1, Path2
        type(OS_type)   :: OS

        call OS%query()
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        Path1 = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path1%Err%occurred
        if (.not. assertion) return

        Path2 = Path_type(inputPath="./temp\ Folder/\{inside\}\/", OS = OS)
        assertion = .not. Path2%Err%occurred
        if (.not. assertion) return

        assertion = Path1%modified == Path2%modified

    end function test_query_3

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
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original  :", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%modified  :", Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%shellSlash:", Path%shellSlash
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir       :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name      :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext       :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_winify_1() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = "./temp\ Folder/\{inside\}/"
        Path%modified = winify(Path%original)

        assertion = Path%modified == '".\temp Folder\{inside}\"'
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") Path%original
            write(Test%outputUnit,"(*(g0,:,' '))") Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_winify_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the routine can successfully remove multiple backslashes from the path to convert them all to a single slash.
    function test_winify_2() result(assertion)

        implicit none
        logical                     :: assertion
        character(:), allocatable   :: modified
        character(:), allocatable   :: original
        original = "\\\"
        modified = winify(original)
        assertion = modified == "\"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") original
            write(Test%outputUnit,"(*(g0,:,' '))") modified
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_winify_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether the routine can successfully convert a single forward-slash path to a backslash.
    function test_winify_3() result(assertion)

        implicit none
        logical                     :: assertion
        character(:), allocatable   :: modified
        character(:), allocatable   :: original
        original = "/"
        modified = winify(original)
        assertion = modified == "\"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") original
            write(Test%outputUnit,"(*(g0,:,' '))") modified
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_winify_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_linify_1() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        Path = Path_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        Path%original = '".\temp Folder\{inside}\"'
        Path%modified = linify(Path%original)

        assertion = Path%modified == "./temp\ Folder/\{inside\}/"

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original     :", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%modified     :", Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%shellSlash   :", Path%shellSlash
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir          :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name         :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext          :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

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
        Path%modified = linify(Path%original)

        assertion = Path%modified == "./temp\ Folder/\{inside\}/-"

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original     :", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%modified     :", Path%modified
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%shellSlash   :", Path%shellSlash
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir          :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name         :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext          :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

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
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirNameExt_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirNameExt_2

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
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirFullName_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When the filename is all file name without extension (any dots), `getDirFullName()` must return
    !> the full file name with empty directory.
    function test_getDirFullName_2() result(assertion)

        implicit none
        logical         :: assertion
        type(Path_type) :: Path

        assertion = .true.

        Path = Path_type(inputPath="ParaMonte")
        assertion = .not. Path%Err%occurred
        if (.not. assertion) return

        call Path%getDirFullName(Path%original,"\",Path%dir,Path%name)

        assertion = assertion .and. Path%dir == ""
        assertion = assertion .and. Path%name == Path%original

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%original :", Path%original
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%dir      :", Path%dir
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name     :", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirFullName_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When the filename is all file extension, `getNameExt()` must return
    !> an empty file name and an extension equivalent to full file name.
    function test_getNameExt_1() result(assertion)

        implicit none
        logical                     :: assertion
        character(:), allocatable   :: name, ext, filename

        assertion = .true.

        filename = ".ParaMonte"
        call getNameExt(filename, name, ext)

        assertion = assertion .and. name == ""
        assertion = assertion .and. ext == filename

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "filename  :", filename
            write(Test%outputUnit,"(*(g0,:,' '))")   "name      :", name
            write(Test%outputUnit,"(*(g0,:,' '))")   "ext       :", ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getNameExt_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getNameExt_2() result(assertion)

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
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(Test%outputUnit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getNameExt_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When the filename is all file extension, `getNameExt()` must return
    !> an empty file name and an extension equivalent to full file name.
    function test_modify_1() result(assertion)

        use System_mod, only: OS_type
        use Err_mod, only: Err_type
        implicit none
        logical                     :: assertion
        character(:), allocatable   :: original, modified, modified_ref
        type(Err_type)              :: Err
        type(OS_type)               :: OS

        original = ".\ParaMonte\dir1 \"
        call modify(original,modified,Err)
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call OS%query()
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (OS%Shell%isUnix) then
            modified_ref = "./ParaMonte/dir1\ /"
            assertion = assertion .and. modified == modified_ref
#if defined OS_IS_WINDOWS
        else
            modified_ref = modified
            assertion = assertion .and. modified == modified
#endif
        end if

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))")   "modified_ref  :", '"'//modified_ref//'"'
            write(Test%outputUnit,"(*(g0,:,' '))")   "modified      :", '"'//modified//'"'
            write(Test%outputUnit,"(*(g0,:,' '))")   "original      :", '"'//original//'"'
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_modify_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether all processors are capable of generating directories.
    function test_mkdir_1() result(assertion)

        use Constants_mod, only: RK
        use System_mod, only: RandomFileName_type, OS_type
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN
        type(OS_type)               :: OS

        RFN = RandomFileName_type(key="test_mkdir_1")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        call OS%query()
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        RFN%Err = mkdir(RFN%path, OS%Shell%isUnix, .true.)
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

    end function test_mkdir_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether all processors are capable of generating directories.
    function test_mkdir_2() result(assertion)

        use Constants_mod, only: RK
        use System_mod, only: RandomFileName_type, OS_type
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN
        type(OS_type)               :: OS

        RFN = RandomFileName_type(key="test_mkdir_2")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        call OS%query()
        assertion = .not. OS%Err%occurred
        if (.not. assertion) return

        RFN%Err = mkdir(RFN%path, OS%Shell%isUnix, .false.)
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

    end function test_mkdir_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test whether all processors are capable of generating directories, without the input optional arguments
    function test_mkdir_3() result(assertion)

        use System_mod, only: RandomFileName_type
        use Constants_mod, only: RK
        implicit none
        logical                     :: assertion
        type(RandomFileName_type)   :: RFN

        RFN = RandomFileName_type(key="test_mkdir_3")
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

        RFN%Err = mkdir(RFN%path)
        assertion = .not. RFN%Err%occurred
        if (.not. assertion) return

    end function test_mkdir_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Path_mod ! LCOV_EXCL_LINE