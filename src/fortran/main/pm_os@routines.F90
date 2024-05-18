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
!>  This file contains procedure implementations of [pm_sysShell](@ref pm_sysShell).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_os) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_sysShell, only: isFailedGetEnvVar
    use pm_sysShell, only: isFailedGetOutput
    use pm_strASCII, only: setStrLower
    use pm_io, only: setContentsFrom
    use pm_io, only: LEN_IOMSG

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isLinux
#if     LINUX_ENABLED
        itis = .true._LK
#else
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isLinux(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isLinux(): "//trim(errmsg) ! LCOV_EXCL_LINE
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isLinuxFailed
#if     LINUX_ENABLED
        failed = .false._LK
        itis = .true._LK
#else
        character(LEN_IOMSG, SK) :: errmsg
        itis = isLinux(failed, errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isLinuxFailedMsg
#if     LINUX_ENABLED
        failed = .false._LK
        itis = .true._LK
#else
        character(:, SK), allocatable :: output
        integer(IK) :: exitstat
        itis = .false._LK
        failed = isFailedGetOutput(SK_"uname - o", output, errmsg = errmsg, exitstat = exitstat)
        if (.not. failed) then
            call setStrLower(output)
            itis = 0 < index(output, "linux")
        end if
#endif
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isDarwin
#if     DARWIN_ENABLED
        itis = .true._LK
#else
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isDarwin(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isDarwin(): "//trim(errmsg) ! LCOV_EXCL_LINE
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isDarwinFailed
#if     DARWIN_ENABLED
        failed = .false._LK
        itis = .true._LK
#else
        character(LEN_IOMSG, SK) :: errmsg
        itis = isDarwin(failed, errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isDarwinFailedMsg
#if     DARWIN_ENABLED
        failed = .false._LK
        itis = .true._LK
#else
        character(:, SK), allocatable :: output
        integer(IK) :: exitstat
        itis = .false._LK
        failed = isFailedGetOutput(SK_"uname - o", output, errmsg = errmsg, exitstat = exitstat)
        if (.not. failed) then
            call setStrLower(output)
            itis = 0 < index(output, "darwin") .or. 0 < index(output, "macos")
        end if
#endif
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isWindows
#if     WINDOWS_ENABLED
        itis = .true._LK
#else
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        itis = isWindows(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isWindows(): "//trim(errmsg) ! LCOV_EXCL_LINE
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isWindowsFailed
#if     WINDOWS_ENABLED
        failed = .false._LK
        itis = .true._LK
#else
        character(LEN_IOMSG, SK) :: errmsg
        itis = isWindows(failed, errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isWindowsFailedMsg
#if     WINDOWS_ENABLED
        failed = .false._LK
        itis = .true._LK
#else
        character(:, SK), allocatable :: value
        itis = .false._LK
        failed = isFailedGetEnvVar(SK_"OS", value, errmsg, length = 10_IK)
        if (.not. failed) then
            call setStrLower(value)
            itis = 0 < index(value, "windows")
        end if
#endif
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines