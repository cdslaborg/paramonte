!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [pm_sysPath](@ref pm_sysPath).
!>  \author Amir Shahmoradi

module test_pm_sysPath

    use pm_sysPath
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)
        call test%run(test_isDir_1, SK_"test_isDir_1")
        call test%run(test_query_1, SK_"test_query_1")
        call test%run(test_query_2, SK_"test_query_2")
        call test%run(test_query_3, SK_"test_query_3")
        call test%run(test_mkdir_1, SK_"test_mkdir_1")
        call test%run(test_mkdir_2, SK_"test_mkdir_2")
        call test%run(test_mkdir_3, SK_"test_mkdir_3")
        call test%run(test_pmodify_1, SK_"test_pmodify_1")
        call test%run(test_getPathWindows_1, SK_"test_getPathWindows_1")
        call test%run(test_getPathWindows_2, SK_"test_getPathWindows_2")
        call test%run(test_getPathWindows_3, SK_"test_getPathWindows_3")
        call test%run(test_getPathPosix_1, SK_"test_getPathPosix_1")
        call test%run(test_getPathPosix_2, SK_"test_getPathPosix_2")
        call test%run(test_getNameExt_1, SK_"test_getNameExt_1")
        call test%run(test_getNameExt_2, SK_"test_getNameExt_2")
        call test%run(test_constructPath, SK_"test_constructPath")
        call test%run(test_getDirNameExt_1, SK_"test_getDirNameExt_1")
        call test%run(test_getDirNameExt_2, SK_"test_getDirNameExt_2")
        call test%run(test_getDirFullName_1, SK_"test_getDirFullName_1")
        call test%run(test_getDirFullName_2, SK_"test_getDirFullName_2")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test if `isDir()` can successfully detect an existing directory.
    function test_isDir_1() result(assertion)
        use pm_kind, only: RK
        implicit none
        logical(LK)     :: assertion
        assertion = isDir("../")
    end function test_isDir_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> When the original path is not allocated, `guery()` must return an error message.
    function test_query_1() result(assertion)

        !use pm_val2str, only: getStr
        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        Path = PathPart_type(inputPath="")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        deallocate(Path%original)
        call Path%query()
        assertion = Path%err%occurred

    end function test_query_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> When the original path is allocated but is empty, `guery()` must return an error message.
    function test_query_2() result(assertion)

        !use pm_val2str, only: getStr
        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        Path = PathPart_type(inputPath="")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = ""
        call Path%query()
        assertion = Path%err%occurred

    end function test_query_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> When the optional OS is provided, the results must be the same as when it is not provided.
    function test_query_3() result(assertion)

        use pm_sysShell, only: OS_type
        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path1, Path2
        type(OS_type)   :: OS

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

        Path1 = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path1%err%occurred
        call test%assert(assertion)

        Path2 = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/", OS = OS)
        assertion = .not. Path2%err%occurred
        call test%assert(assertion)

        assertion = Path1%modified == Path2%modified

    end function test_query_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPath() result(assertion)

        !use pm_val2str, only: getStr
        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        assertion = Path%original == "./temp\ Folder/\{inside\}\/"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%original  :", Path%original
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%modified  :", Path%modified
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%sep:", Path%sep
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir       :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name      :", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%ext       :", Path%ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPathWindows_1() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = "./temp\ Folder/\{inside\}/"
        Path%modified = getPathWindows(Path%original)

        assertion = Path%modified == '".\temp Folder\{inside}\"'
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") Path%original
            write(test%disp%unit,"(*(g0,:,' '))") Path%modified
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPathWindows_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether the routine can successfully remove multiple backslashes from the path to convert them all to a single slash.
    function test_getPathWindows_2() result(assertion)

        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: modified
        character(:, SK), allocatable   :: original
        original = "\\\"
        modified = getPathWindows(original)
        assertion = modified == "\"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") original
            write(test%disp%unit,"(*(g0,:,' '))") modified
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPathWindows_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether the routine can successfully convert a single forward-slash path to a backslash.
    function test_getPathWindows_3() result(assertion)

        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: modified
        character(:, SK), allocatable   :: original
        original = "/"
        modified = getPathWindows(original)
        assertion = modified == "\"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))") original
            write(test%disp%unit,"(*(g0,:,' '))") modified
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPathWindows_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPathPosix_1() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = '".\temp Folder\{inside}\"'
        Path%modified = getPathPosix(Path%original)

        assertion = Path%modified == "./temp\ Folder/\{inside\}/"

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%original     :", Path%original
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%modified     :", Path%modified
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%sep   :", Path%sep
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir          :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name         :", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%ext          :", Path%ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPathPosix_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPathPosix_2() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = '".\temp Folder\{inside}\-"'
        Path%modified = getPathPosix(Path%original)

        assertion = Path%modified == "./temp\ Folder/\{inside\}/-"

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%original     :", Path%original
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%modified     :", Path%modified
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%sep   :", Path%sep
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir          :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name         :", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%ext          :", Path%ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPathPosix_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDirNameExt_1() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        assertion = .true._LK

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = ".\temp Folder\{inside}\"
        call Path%getDirNameExt(Path%original,"\",Path%dir,Path%name,Path%ext)

        assertion = assertion .and. Path%dir == ".\temp Folder\{inside}\"
        assertion = assertion .and. Path%name == ""
        assertion = assertion .and. Path%ext == ""

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirNameExt_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDirNameExt_2() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        assertion = .true._LK

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = ".\temp Folder\{inside}\Temp.tXt"
        call Path%getDirNameExt(Path%original,"\",Path%dir,Path%name,Path%ext)

        assertion = assertion .and. Path%dir == ".\temp Folder\{inside}\"
        assertion = assertion .and. Path%name == "Temp"
        assertion = assertion .and. Path%ext == ".tXt"

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirNameExt_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getDirFullName_1() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        assertion = .true._LK

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = ".\temp Folder\{inside}\-"
        call Path%getDirFullName(Path%original,"\",Path%dir,Path%name)

        assertion = assertion .and. Path%dir == ".\temp Folder\{inside}\"
        assertion = assertion .and. Path%name == "-"

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirFullName_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> When the filename is all file name without extension (any dots), `getDirFullName()` must return
    !> the full file name with empty directory.
    function test_getDirFullName_2() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        assertion = .true._LK

        Path = PathPart_type(inputPath="ParaMonte")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        call Path%getDirFullName(Path%original,"\",Path%dir,Path%name)

        assertion = assertion .and. Path%dir == ""
        assertion = assertion .and. Path%name == Path%original

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%original :", Path%original
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%dir      :", Path%dir
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name     :", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getDirFullName_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> When the filename is all file extension, `getNameExt()` must return
    !> an empty file name and an extension equivalent to full file name.
    function test_getNameExt_1() result(assertion)

        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: name, ext, filename

        assertion = .true._LK

        filename = SK_".ParaMonte"
        call getNameExt(filename, name, ext)

        assertion = assertion .and. name == ""
        assertion = assertion .and. ext == filename

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "filename  :", filename
            write(test%disp%unit,"(*(g0,:,' '))")   "name      :", name
            write(test%disp%unit,"(*(g0,:,' '))")   "ext       :", ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getNameExt_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getNameExt_2() result(assertion)

        implicit none
        logical(LK)     :: assertion
        type(PathPart_type) :: Path

        assertion = .true._LK

        Path = PathPart_type(inputPath="./temp\ Folder/\{inside\}\/")
        assertion = .not. Path%err%occurred
        call test%assert(assertion)

        Path%original = "-"
        call Path%getNameExt(Path%original,Path%name,Path%ext)

        assertion = assertion .and. Path%name == "-"
        assertion = assertion .and. Path%ext == ""

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%name:", Path%name
            write(test%disp%unit,"(*(g0,:,' '))")   "Path%ext :", Path%ext
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getNameExt_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> When the filename is all file extension, `getNameExt()` must return
    !> an empty file name and an extension equivalent to full file name.
    function test_pmodify_1() result(assertion)

        use pm_sysShell, only: OS_type
        use pm_err, only: err_type
        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: original, modified, modified_ref
        type(err_type)              :: Err
        type(OS_type)               :: OS

        original = ".\ParaMonte\dir1 \"
        call modify(original,modified,Err)
        assertion = .not. OS%err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        call OS%query()
        assertion = .not. OS%err%occurred
        if (.not. assertion) return ! LCOV_EXCL_LINE

        if (OS%shell%is%posix) then
            modified_ref = "./ParaMonte/dir1\ /"
            assertion = assertion .and. modified == modified_ref
#if WINDOWS_ENABLED
        else
            modified_ref = modified
            assertion = assertion .and. modified == modified
#endif
        end if

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,' '))")
            write(test%disp%unit,"(*(g0,:,' '))")   "modified_ref  :", '"'//modified_ref//'"'
            write(test%disp%unit,"(*(g0,:,' '))")   "modified      :", '"'//modified//'"'
            write(test%disp%unit,"(*(g0,:,' '))")   "original      :", '"'//original//'"'
            write(test%disp%unit,"(*(g0,:,' '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_pmodify_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether all processors are capable of generating directories.
    function test_mkdir_1() result(assertion)

        use pm_kind, only: RK
        use pm_sysShell, only: OS_type
        use pm_sysPath, only: getPathNew
        logical(LK) :: assertion
        character(:, SK), allocatable   :: newpath
        type(OS_type)                   :: OS

        newpath = getPathNew(prefix = SK_"test_mkdir_1", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

        assertion = .not. isFailedMakeDir(newpath)
        call test%assert(assertion)

    end function test_mkdir_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether all processors are capable of generating directories.
    function test_mkdir_2() result(assertion)

        use pm_kind, only: RK
        use pm_sysShell, only: OS_type
        use pm_sysPath, only: getPathNew
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: newpath
        type(OS_type)                   :: OS

        newpath = getPathNew(prefix = SK_"test_mkdir_2", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

        assertion = .not. isFailedMakeDir(newpath)
        call test%assert(assertion)

    end function test_mkdir_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether all processors are capable of generating directories, without the input optional arguments
    function test_mkdir_3() result(assertion)

        use pm_err, only: err_type
        use pm_kind, only: RK
        use pm_sysPath, only: getPathNew
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: newpath
        type(err_type)                  :: Err

        newpath = getPathNew(prefix = SK_"test_mkdir_3", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        assertion = .not. isFailedMakeDir(newpath)
        call test%assert(assertion)

    end function test_mkdir_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_sysPath ! LCOV_EXCL_LINE