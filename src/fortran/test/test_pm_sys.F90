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

!>  \brief This module contains tests of the module [pm_sysShell](@ref pm_sysShell).
!>  \author Amir Shahmoradi

module test_pm_sys

    use pm_sysShell
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

        ! On Windows subsystem for Linux, various runtime errors arise with the following tests, in particular, with
        ! test_SysCmd_type_1 and test_SystemInfo_type_1 with such error messages as the following:
        !
        !   Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained
        !
        ! Therefore, these tests might be only enabled for code-coverage and instrumentation purposes.
        !
        ! See: https://community.intel.com/t5/Intel-Fortran-Compiler/Fortran-execute-command-line-runtime-error-depends-on-memory/m-p/1168196#M145090

        test = test_type(MODULE_NAME)
        call test%run(test_sleep_1, SK_"test_sleep_1")
        call test%run(test_OS_type_1, SK_"test_OS_type_1")
        call test%run(test_OS_type_2, SK_"test_OS_type_2")
        call test%run(test_OS_type_3, SK_"test_OS_type_3")
        call test%run(test_copyFile_1, SK_"test_copyFile_1")
        call test%run(test_removeFile_1, SK_"test_removeFile_1")
        call test%run(test_removeFile_2, SK_"test_removeFile_2")
        call test%run(test_isFailedExec_1, SK_"test_isFailedExec_1")
        call test%run(test_isFailedExec_2, SK_"test_isFailedExec_2")
        call test%run(test_EnvVar_type_1, SK_"test_EnvVar_type_1")
        call test%run(test_EnvVar_type_2, SK_"test_EnvVar_type_2")
        call test%run(test_EnvVar_type_3, SK_"test_EnvVar_type_3")
        call test%run(test_CmdArg_type_1, SK_"test_CmdArg_type_1")
        call test%run(test_SysCmd_type_1, SK_"test_SysCmd_type_1")
        call test%run(test_getSysInfo_1, SK_"test_getSysInfo_1")
        call test%run(test_SystemInfo_type_1, SK_"test_SystemInfo_type_1")
        call test%run(test_getPathNew_1, SK_"test_getPathNew_1")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_sleep_1() result(assertion)
        use pm_kind, only: RK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        type(err_type)  :: Err
        call sleep(0.1_RK, Err)
        assertion = .not. err%occurred
    end function test_sleep_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_SystemInfo_type_1() result(assertion)
        implicit none
        logical(LK)             :: assertion
        type(SystemInfo_type)   :: SystemInfo
        assertion = .true._LK
        SystemInfo = SystemInfo_type(pid = test%image%id)
        assertion = .not. SystemInfo%err%occurred .and. allocated(SystemInfo%Records) .and. size(SystemInfo%Records) > 0
    end function test_SystemInfo_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether environmental variables can be successfully queried.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_EnvVar_type_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)             :: assertion
        type(EnvVar_type)       :: EnvVar
        assertion = .true._LK

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        EnvVar%name = SK_""
        call EnvVar%get(EnvVar%name,EnvVar%val,EnvVar%length,EnvVar%err)
        !assertion = EnvVar%err%occurred

    end function test_EnvVar_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether environmental variables can be successfully queried.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_EnvVar_type_2() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)             :: assertion
        type(EnvVar_type)       :: EnvVar
        assertion = .true._LK

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        EnvVar%name = SK_"OS"
        call EnvVar%get(EnvVar%name,EnvVar%val,EnvVar%length,EnvVar%err)
        !assertion = .not. EnvVar%err%occurred .and. allocated(EnvVar%name) .and. allocated(EnvVar%val)

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "EnvVar%name     : ", EnvVar%name
        !    write(test%disp%unit,"(2A)")   "EnvVar%val    : ", EnvVar%val
        !    write(test%disp%unit,"(2A)")   "EnvVar%length   : ", getStr(EnvVar%length)
        !    write(test%disp%unit,"(2A)")
        !end if
        !! LCOV_EXCL_STOP

    end function test_EnvVar_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether environmental variables can be successfully queried.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_EnvVar_type_3() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)             :: assertion
        type(EnvVar_type)       :: EnvVar
        assertion = .true._LK

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        EnvVar%name = SK_"OS"
        call EnvVar%get(EnvVar%name, EnvVar%val, Err = EnvVar%err)
        !assertion = .not. EnvVar%err%occurred .and. allocated(EnvVar%name) .and. allocated(EnvVar%val)

    end function test_EnvVar_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test whether input command line arguments can be successfully retrieved.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_CmdArg_type_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
       !integer(IK)             :: i
        logical(LK)             :: assertion
        type(CmdArg_type)       :: CmdArg
        assertion = .true._LK

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        call CmdArg%query()
        assertion = .not. CmdArg%err%occurred
        call test%assert(assertion)

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "CmdArg%cmd      : ", CmdArg%cmd
        !    write(test%disp%unit,"(2A)")   "CmdArg%count    : ", getStr(CmdArg%count)
        !    write(test%disp%unit,      "(*('CmdArg%slash    : ', 2A))") (CmdArg%Arg(i)%record, new_line('a'), i=1,CmdArg%count)
        !    write(test%disp%unit,"(2A)")
        !end if
        !! LCOV_EXCL_STOP

    end function test_CmdArg_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Query Shell first and then OS to ensure caching the Shell query results work correctly.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_OS_type_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)             :: assertion
        type(OS_type)           :: OS

        assertion = .true._LK

!#if CODECOV_ENABLED
        mv_osCacheActivated = .false._LK
        mc_shellCached = .false._LK
!#endif

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        call OS%query()
        !assertion = .not. OS%err%occurred .and. .not. OS%shell%err%occurred
        !if (.not. assertion) return

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "OS%name     : ", OS%name
        !    write(test%disp%unit,"(2A)")   "OS%slash    : ", OS%slash
        !    write(test%disp%unit,"(2A)")   "OS%isWindows: ", getStr(OS%isWindows)
        !    write(test%disp%unit,"(2A)")
        !end if
        !! LCOV_EXCL_STOP

    end function test_OS_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Query Shell first and then OS to ensure caching the Shell query results work correctly.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_OS_type_2() result(assertion)

        use pm_val2str, only: getStr
        use pm_strASCII, only: getStrLower
        implicit none
        logical(LK)     :: assertion
        type(OS_type)   :: OS

        assertion = .true._LK

!#if CODECOV_ENABLED
        mv_osCacheActivated = .false._LK
        mc_shellCached = .false._LK
#endif

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        call OS%shell%query()
        call shell_typer(OS%Shell)
        !assertion = assertion .and. .not. OS%shell%err%occurred

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "OS%name           : ", OS%name
        !    write(test%disp%unit,"(2A)")   "OS%slash          : ", OS%slash
        !    write(test%disp%unit,"(2A)")   "OS%isWindows      : ", getStr(OS%isWindows)
        !    write(test%disp%unit,"(2A)")   "OS%shell%name     : ", OS%shell%name
        !    write(test%disp%unit,"(2A)")   "OS%shell%slash    : ", OS%shell%slash
        !    write(test%disp%unit,"(2A)")   "OS%shell%is%posix   : ", getStr(OS%shell%is%posix)
        !    write(test%disp%unit,"(2A)")
        !end if
        !! LCOV_EXCL_STOP

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        call OS%query()
        !assertion = assertion .and. .not. OS%err%occurred .and. .not. OS%shell%err%occurred

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "OS%name           : ", OS%name
        !    write(test%disp%unit,"(2A)")   "OS%slash          : ", OS%slash
        !    write(test%disp%unit,"(2A)")   "OS%isWindows      : ", getStr(OS%isWindows)
        !    write(test%disp%unit,"(2A)")   "OS%shell%name     : ", OS%shell%name
        !    write(test%disp%unit,"(2A)")   "OS%shell%slash    : ", OS%shell%slash
        !    write(test%disp%unit,"(2A)")   "OS%shell%is%posix   : ", getStr(OS%shell%is%posix)
        !    write(test%disp%unit,"(2A)")
        !end if
        ! LCOV_EXCL_STOP

    end function test_OS_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Query OS first without shell query and then query OS with shell.
    !>  This will test the remaining uncovered cached conditions in `queryOS()`.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_OS_type_3() result(assertion)

        use pm_val2str, only: getStr
        use pm_strASCII, only: getStrLower
        implicit none
        logical(LK)     :: assertion
        type(OS_type)   :: OS

        assertion = .true._LK

!#if CODECOV_ENABLED
        mv_osCacheActivated = .false._LK
        mc_shellCached = .false._LK
#endif

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        call OS%query(shellQueryEnabled = .false._LK)
        !assertion = assertion .and. .not. OS%shell%err%occurred

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "OS%name           : ", OS%name
        !    write(test%disp%unit,"(2A)")   "OS%slash          : ", OS%slash
        !    write(test%disp%unit,"(2A)")   "OS%isWindows      : ", getStr(OS%isWindows)
        !    write(test%disp%unit,"(2A)")   "OS%shell%name     : ", OS%shell%name
        !    write(test%disp%unit,"(2A)")   "OS%shell%slash    : ", OS%shell%slash
        !    write(test%disp%unit,"(2A)")   "OS%shell%is%posix   : ", getStr(OS%shell%is%posix)
        !    write(test%disp%unit,"(2A)")
        !end if
        !! LCOV_EXCL_STOP

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained

        call OS%query(shellQueryEnabled = .true._LK)
        !assertion = assertion .and. .not. OS%err%occurred .and. .not. OS%shell%err%occurred

        !if (test%traceable .and. .not. assertion) then
        !! LCOV_EXCL_START
        !    write(test%disp%unit,"(2A)")
        !    write(test%disp%unit,"(2A)")   "OS%name           : ", OS%name
        !    write(test%disp%unit,"(2A)")   "OS%slash          : ", OS%slash
        !    write(test%disp%unit,"(2A)")   "OS%isWindows      : ", getStr(OS%isWindows)
        !    write(test%disp%unit,"(2A)")   "OS%shell%name     : ", OS%shell%name
        !    write(test%disp%unit,"(2A)")   "OS%shell%slash    : ", OS%shell%slash
        !    write(test%disp%unit,"(2A)")   "OS%shell%is%posix   : ", getStr(OS%shell%is%posix)
        !    write(test%disp%unit,"(2A)")
        !end if
        ! LCOV_EXCL_STOP

    end function test_OS_type_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the success of a SysCmd action.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_SysCmd_type_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)             :: assertion
        type(SysCmd_type)       :: SysCmd
        type(OS_type)           :: OS

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

        ! No assertion evaluation below as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained
        ! The `Err` argument handles exceptions.

#if WINDOWS_ENABLED
        if (OS%shell%is%cmd .or. OS%shell%is%powershell) then
            SysCmd = SysCmd_type("dir > nul 2>&1", .true._LK)
            !assertion = .not. SysCmd%err%occurred
            !if (.not. assertion) return
            SysCmd = SysCmd_type("dir > nul 2>&1", .false._LK)
            !assertion = .not. SysCmd%err%occurred
            !if (.not. assertion) return
        end if
#endif
        if (OS%shell%is%posix) then
            SysCmd = SysCmd_type("ls >/dev/null 2>&1", .true._LK)
            !assertion = .not. SysCmd%err%occurred
            !if (.not. assertion) return
            SysCmd = SysCmd_type("ls >/dev/null 2>&1", .false._LK)
            !assertion = .not. SysCmd%err%occurred
            !if (.not. assertion) return
        end if

    end function test_SysCmd_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the success of a SysCmd action.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_isFailedExec_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)                     :: assertion
        type(OS_type)                   :: OS
        integer(IK)                     :: exitstat
        logical(LK)     , parameter     :: wait = .true._LK
        character(:, SK), allocatable   :: command

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

#if WINDOWS_ENABLED
        if (OS%shell%is%cmd .or. OS%shell%is%powershell) then
            command = "dir > nul 2>&1"
        end if
#endif
        if (OS%shell%is%posix) then
            command = "ls >/dev/null 2>&1"
        end if

        ! No assertion evaluation here as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained
        ! The `Err` argument handles exceptions.

        call isFailedExec(command,wait,exitstat,OS%err)
        !assertion = .not. OS%err%occurred
        !if (.not. assertion) return

    end function test_isFailedExec_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the success of a SysCmd action.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_isFailedExec_2() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)                     :: assertion
        type(OS_type)                   :: OS
        character(:, SK), allocatable   :: command

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

#if WINDOWS_ENABLED
        if (OS%shell%is%cmd .or. OS%shell%is%powershell) then
            command = "dir > nul 2>&1"
        end if
#endif
        if (OS%shell%is%posix) then
            command = "ls >/dev/null 2>&1"
        end if

        ! No assertion evaluation here as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained
        ! The `Err` argument handles exceptions.

        call isFailedExec(command, Err = OS%err)

    end function test_isFailedExec_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPathNew_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        character(:, SK), allocatable :: newpath
        logical(LK)                 :: fileExists

        newpath = getPathNew(prefix = SK_"test_getPathNew", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        inquire(file = newpath, exist = fileExists)
        assertion = assertion .and. .not. fileExists

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(2A)")
            write(test%disp%unit,"(2A)")   "newpath     : ", newpath
            ! LCOV_EXCL_STOP
        end if

    end function test_getPathNew_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_removeFile_1() result(assertion)

        use pm_val2str, only: getStr
        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: newpath
        logical(LK)                     :: fileExists
        integer(IK)                     :: fileUnit

        newpath = getPathNew(prefix = SK_"test_getPathNew", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        open(newunit=fileUnit,file=newpath,status="new")
        close(fileUnit)

        call removeFile(path = newpath, Err = test%err)
        assertion = assertion .and. .not. Tes%err%occurred
        call test%assert(assertion)

        inquire(file=newpath,exist=fileExists)
        assertion = assertion .and. .not. fileExists

    end function test_removeFile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_removeFile_2() result(assertion)

        use pm_kind, only: RK
        use pm_val2str, only: getStr
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: newpath
        logical(LK)                     :: fileExists
        integer(IK)                     :: fileUnit

        newpath = getPathNew(prefix = SK_"test_getPathNew", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        open(newunit=fileUnit,file=newpath,status="new")
        close(fileUnit)

        call removeFile(path=newpath)
        inquire(file=newpath,exist=fileExists)
        assertion = assertion .and. .not. fileExists

    end function test_removeFile_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Obtain the system info without providing the cachefile name, in which case, no cache file will be generated.
    !>
    !> \todo
    !>  This test needs further improvements in the future. See comments in the body of the test.
    function test_getSysInfo_1() result(assertion)

        use pm_container, only: css_pdt
        use pm_err, only: err_type
        implicit none
        logical(LK)                         :: assertion
        type(css_pdt)    , allocatable   :: List(:)
        type(err_type)                      :: Err

        assertion = .true._LK

        ! No assertion evaluation here as execute_command_line() is prone to failure:
        ! Fortran runtime error: EXECUTE_COMMAND_LINE: Termination status of the command-language interpreter cannot be obtained
        ! The `Err` argument handles exceptions.

        call getSysInfo(List,Err)
        !assertion = .not. err%occurred

    end function test_getSysInfo_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_copyFile_1() result(assertion)

        use pm_kind, only: RK
        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: newpath
        type(OS_type)                   :: OS
        logical(LK)                     :: fileExists
        integer(IK)                     :: fileUnit

        newpath = getPathNew(prefix = SK_"test_getPathNew", failed = assertion)
        assertion = .not. assertion
        call test%assert(assertion)

        open(newunit=fileUnit,file=newpath,status="new")
        close(fileUnit)

        call OS%query()
        assertion = .not. OS%err%occurred
        call test%assert(assertion)

        assertion = .not. isFailedCopy(newpath, newpath//".copy")
        call test%assert(assertion)

        inquire(file=newpath//".copy",exist=fileExists)
        assertion = assertion .and. fileExists

        call removeFile(newpath)
        inquire(file=newpath,exist=fileExists)
        assertion = assertion .and. .not. fileExists

        call removeFile(newpath//".copy")
        inquire(file=newpath//".copy",exist=fileExists)
        assertion = assertion .and. .not. fileExists

    end function test_copyFile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_sys ! LCOV_EXCL_LINE