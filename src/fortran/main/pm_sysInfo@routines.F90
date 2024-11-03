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

!>  \brief
!>  This file contains procedure implementations of [pm_sysInfo](@ref pm_sysInfo).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
!   The Intel FPP defines `linux=1` which conflicts with the `linux` component of `kernelis_type` below.
#if     __INTEL_COMPILER
#undef  linux
#endif
!>  \endcond excluded

submodule (pm_sysInfo) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_io, only: LEN_IOMSG
    use pm_strASCII, only: getStrLower
    use pm_sysShell, only: isFailedGetEnvVar
    use pm_sysShell, only: isFailedGetOutput

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `private` module object of type [kernel_type](@ref pm_sysInfo::kernel_type).
    !>
    !>  \details
    !>  This object is set only once throughout the life of the program to avoid costly redundant construction of OS objects.<br>
    !>  It is set by and exclusively used within the routines of this submodule and nowhere else.<br>
    !>  The allocation status of the object is used as an indicator of its initialization.<br>
    !>
    !>  \final{kernel_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(kernel_type), allocatable :: mc_kernel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function isFailedInitOS(errmsg) result(failed)
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(inout), optional :: errmsg
        logical(LK) :: failed
        failed = .not. allocated(mc_kernel)
        if (failed) then
            if (present(errmsg)) then
                mc_kernel = kernel_type(failed, errmsg)
            else
                mc_kernel = kernel_type(failed)
            end if
            if (failed) deallocate(mc_kernel)
        end if
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure kernelis_typer
        character(255, SK) :: errmsg
        logical(LK) :: failed
        errmsg = SK_""
        kernelis = kernelis_type(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@kernelis_typer(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure kernelis_typerFailed
        character(255, SK) :: errmsg
        errmsg = SK_""
        kernelis = kernelis_type(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure kernelis_typerFailedMsg
        if (allocated(mc_kernel)) then
            kernelis = mc_kernel%is
            failed = .false._LK
        else
            mc_kernel = kernel_type(failed, errmsg)
            if (allocated(mc_kernel)) then
                kernelis = mc_kernel%is
            else
                errmsg = MODULE_NAME//SK_"@kernelis_typerFailedMsg(): "//trim(errmsg) ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            end if
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure kernel_typer
        character(255, SK) :: errmsg
        logical(LK) :: failed
        errmsg = SK_""
        kernel = kernel_type(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@kernel_typer(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure kernel_typerFailed
        character(31, SK) :: errmsg
        errmsg = SK_""
        kernel = kernel_type(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure kernel_typerFailedMsg
        use pm_kind, only: SKG => SK
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@kernel_typer()"
        if (allocated(mc_kernel)) then
            failed = .false._LK
            kernel = mc_kernel
        else
#if         MSYS_ENABLED
            kernel%name = "MSYS"
            kernel%is%msys = .true._LK
            kernel%is%windows = .true._LK
            failed = .false._LK
#elif       MINGW_ENABLED
            kernel%name = "MinGW"
            kernel%is%mingw = .true._LK
            kernel%is%windows = .true._LK
            failed = .false._LK
#elif       CYGWIN_ENABLED
            kernel%name = "Cygwin"
            kernel%is%cygwin = .true._LK
            kernel%is%windows = .true._LK
            failed = .false._LK
#elif       WINDOWS_ENABLED
            kernel%name = "Windows"
            kernel%is%windows = .true._LK
            failed = .false._LK
#elif       DARWIN_ENABLED
            kernel%name = "Darwin"
            kernel%is%darwin = .true._LK
            failed = .false._LK
#elif       LINUX_ENABLED
            kernel%name = "Linux"
            kernel%is%linux = .true._LK
            failed = .false._LK
#else
            block
                character(:,SKG), allocatable :: nameLower
                ! Lets infer it through the runtime shell. Firs assume Window OS.
                failed = isFailedGetEnvVar(name = SK_"OS", value = kernel%name, errmsg, length = length)
                if (failed) then
                    ! assume unix-like.
                    failed = isFailedGetEnvVar(name = SK_"OSTYPE", value = kernel%name, errmsg, length = length)
                    if (failed) then
                        ! Try one last time with the universal UNIX command `uname`.
                        failed = isFailedGetOutput(SKG_"uname", kernel%name, errmsg)
                        if (failed) then
                            errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg) ! LCOV_EXCL_LINE
                            kernel%name = SKG_"" ! LCOV_EXCL_LINE
                            return ! LCOV_EXCL_LINE
                        end if
                    end if
                end if
                nameLower = getStrLower(kernel%name)
                kernel%is%freebsd = logical(index(nameLower, SKG_"freebsd"), LK)
                kernel%is%darwin = logical(index(nameLower, SKG_"darwin"), LK)
                kernel%is%linux = logical(index(nameLower, SKG_"linux"), LK)
                kernel%is%msys = logical(index(nameLower, SKG_"msys"), LK)
                kernel%is%mingw = logical(index(nameLower, SKG_"mingw"), LK)
                kernel%is%cygwin = logical(index(nameLower, SKG_"cygwin"), LK)
                kernel%is%windows = logical(index(nameLower, SKG_"windows"), LK) .or. kernel%is%msys .or. kernel%is%mingw .or. kernel%is%cygwin
            end block
#endif
            if (.not. allocated(mc_kernel)) allocate(mc_kernel, source = kernel)
            !mc_kernel = kernel
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isKernelWindows
#if     WINDOWS_ENABLED
        itis = .true._LK
#else
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isKernelWindows(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isKernelWindows(): "//trim(errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isKernelWindowsFailed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isKernelWindows(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isKernelWindows_ENABLED 1
    module procedure isKernelWindowsFailedMsg
#include "pm_sysInfo@routines.inc.F90"
    end procedure
#undef isKernelWindows_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isKernelDarwin
#if     DARWIN_ENABLED
        itis = .true._LK
#else
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isKernelDarwin(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isKernelDarwin(): "//trim(errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isKernelDarwinFailed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isKernelDarwin(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isKernelDarwin_ENABLED 1
    module procedure isKernelDarwinFailedMsg
#include "pm_sysInfo@routines.inc.F90"
    end procedure
#undef isKernelDarwin_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isKernelLinux
#if     LINUX_ENABLED
        itis = .true._LK
#else
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isKernelLinux(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isKernelLinux(): "//trim(errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isKernelLinuxFailed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isKernelLinux(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isKernelLinux_ENABLED 1
    module procedure isKernelLinuxFailedMsg
#include "pm_sysInfo@routines.inc.F90"
    end procedure
#undef isKernelLinux_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getSysInfo
        use pm_kind, only: SKG => SK
        character(LEN_IOMSG,SKG) :: errmsg
        logical(LK) :: failed
        errmsg = SKG_""
        sysInfo = getSysInfo(failed, errmsg)
        if (failed) error stop trim(errmsg)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getSysInfoFailed
        use pm_kind, only: SKG => SK
        character(LEN_IOMSG,SKG) :: errmsg
        sysInfo = getSysInfo(failed, errmsg)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getSysInfoFailedMsg

        use pm_kind, only: SKG => SK
        use pm_val2str, only: getStr
        use pm_io, only: setContentsFrom
        use pm_sysShell, only: isFailedExec
        use pm_sysPath, only: isFailedRemove
        use pm_sysPath, only: getPathTemp
        use pm_sysShell, only: shell_type
        use pm_container, only: css_type

        character(*,SKG), parameter     :: NLC = new_line(SKG_"a")
        character(*, SK), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@getSysInfo()"
        character(:,SKG), allocatable   :: contents, stderr, stdout, dumper
        type(css_type)  , allocatable   :: cmd(:)
        type(shell_type)                :: shell
        logical(LK)                     :: done
        integer(IK)                     :: icmd
        integer(IK)                     :: iostat

        sysInfo = SKG_""

        !!!!
        !!!! Infer the shell type.
        !!!!

        shell = shell_type(failed, errmsg)
        if (failed) then
            errmsg = PROCEDURE_NAME//SK_": Failed to create a new file path for system information storage." ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        ! Define the cache file.
        !if (present(cache)) then
        !    cache_def = cache
        !else
        !    ! generate a brand new, non-existing temporary filename.
        !    cache_def = getPathTemp(prefix = SK_".sysInfo", failed = failed)
        !    if (failed) then
        !        if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Failed to create a new file path for system information storage." ! LCOV_EXCL_LINE
        !        return ! LCOV_EXCL_LINE
        !    end if
        !end if

        !!!!
        !!!! Standard error output for error-catching.
        !!!!

       !stderr = SK_" 2> "//cache_def//SK_".stderr"
        stdout = getPathTemp(prefix = SK_".sysInfo", failed = failed)
        stderr = stdout//SK_".stderr"
        dumper = SK_" 1> "//stdout//SK_" 2> "//stderr

        !!!!
        !!!! Define the shell command.
        !!!!

        failed = isFailedInitOS(errmsg)
        if (failed) return ! LCOV_EXCL_LINE

        if (mc_kernel%is%darwin) then
            cmd =   [ css_type(SK_"uname -a"//dumper) &
                    , css_type(SK_"sysctl -a | grep machdep.cpu"//dumper) &
                    , css_type(SK_"system_profiler SPHardwareDataType"//dumper) &
                    ]
        elseif (mc_kernel%is%linux .or. mc_kernel%is%freebsd .or. shell%is%posix) then
            !command = SK_"uname -a > "//cache_def//SK_"; lshw -short >> "//cache_def//SK_"; lscpu >> "//cache_def
            cmd =   [ css_type(SK_"uname -a"//dumper) &
                    , css_type(SK_"lscpu"//dumper//SK_" || cat /proc/cpuinfo"//dumper) &
                    ]
        elseif (mc_kernel%is%windows .and. shell%is%windows) then
            cmd =   [css_type(SK_"systeminfo"//dumper)]
        else
            errmsg = PROCEDURE_NAME//SK_": Unrecognized operating system: "//mc_kernel%name ! LCOV_EXCL_LINE
            failed = .true._LK ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        !!!!
        !!!! Get sysinfo.
        !!!!

        done = .false._LK
        do icmd = 1, size(cmd, 1, IK)
            if (.not. isFailedExec(cmd(icmd)%val, cmdmsg = errmsg)) then
                call setContentsFrom(stdout, contents, iostat, iomsg = errmsg, del = .true._LK)
                if (iostat == 0_IK) then
                    sysInfo = sysInfo//NLC//contents
                    done = .true._LK
                end if
            end if
        end do

        !!!!
        !!!! At least one of the loop cycles must succeed to not fail.
        !!!!

        failed = .not. done
        if (failed) then
            errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg) ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        !   \todo
        !   \warning
        !   On some platforms, such as Windows Subsystem for Linux, the CMD exit status might not
        !   be returned reliably and therefore, cause `isFailedExec()` to return an error.
        !   In such a case, no error for the file copying should be really raised.
        !   If the file already exists upon copy action, no error should be raised.
        !   Note that this method may have some vulnerabilities, for example, when
        !   a file copy is created, but the copy action did not accomplish the
        !   task successfully and the copied file is broken.
        !   This needs a more robust solution in the future.
        !failed = .true._LK ! LCOV_EXCL_LINE
        !return ! LCOV_EXCL_LINE

        !!!!
        !!!! Delete the stderr file.
        !!!!

        failed = isFailedRemove(stderr, errmsg = errmsg)

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
