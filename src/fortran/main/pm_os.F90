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
!>  This module contains procedures and generic interfaces for inferring the processor operating system.
!>
!>  \see
!>  [pm_io](@ref pm_io)<br>
!>  [pm_sysInfo](@ref pm_sysInfo)<br>
!>  [pm_sysPath](@ref pm_sysPath)<br>
!>  [pm_sysShell](@ref pm_sysShell)<br>
!>
!>  \test
!>  [test_pm_os](@ref test_pm_os)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 3:09 AM, Dec 8, 2017, Dell Medical School, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_os

    use pm_kind, only: IK, LK, RK, CK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_os"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime operating system is Linux.
    !>
    !>  \details
    !>  The identification of an operating system as Linux requires either:<br>
    !>  <ol>
    !>      <li>    the ParaMonte library to have been compiled with the FPP macro `LINUX_ENABLED=1`, or,
    !>      <li>    the value returned by the runtime shell command `uname -o` match or contain the keyword `linux`.
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `itis`                  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> the runtime operating system is Linux, otherwise, `.false.`.
    !>
    !>  \interface{isLinux}
    !>  \code{.F90}
    !>
    !>      use pm_os, only: isLinux
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: itis
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      itis = isLinux()
    !>      itis = isLinux(failed)
    !>      itis = isLinux(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isLinux](@ref pm_os::isLinux)<br>
    !>  [isDarwin](@ref pm_os::isDarwin)<br>
    !>  [isWindows](@ref pm_os::isWindows)<br>
    !>
    !>  \example{isLinux}
    !>  \include{lineno} example/pm_os/isLinux/main.F90
    !>  \compilef{isLinux}
    !>  \output{isLinux}
    !>  \include{lineno} example/pm_os/isLinux/main.out.F90
    !>
    !>  \test
    !>  [test_pm_os](@ref test_pm_os)
    !>
    !>  \final{isLinux}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isLinux

    impure module function isLinux() result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLinux
#endif
        logical(LK)                     :: itis
    end function

    impure module function isLinuxFailed(failed) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLinuxFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: itis
    end function

    impure module function isLinuxFailedMsg(failed, errmsg) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isLinuxFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: itis
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime operating system is Darwin (macOS).
    !>
    !>  \details
    !>  The identification of an operating system as Darwin requires either:<br>
    !>  <ol>
    !>      <li>    the ParaMonte library to have been compiled with the FPP macro `DARWIN_ENABLED=1`, or,
    !>      <li>    the value returned by the runtime shell command `uname -o` match or contain the keyword `darwin` or `macos`.
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `itis`                  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> the runtime operating system is Darwin (macOS), otherwise, `.false.`.
    !>
    !>  \interface{isDarwin}
    !>  \code{.F90}
    !>
    !>      use pm_os, only: isDarwin
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: itis
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      itis = isDarwin()
    !>      itis = isDarwin(failed)
    !>      itis = isDarwin(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isLinux](@ref pm_os::isLinux)<br>
    !>  [isDarwin](@ref pm_os::isDarwin)<br>
    !>  [isWindows](@ref pm_os::isWindows)<br>
    !>
    !>  \example{isDarwin}
    !>  \include{lineno} example/pm_os/isDarwin/main.F90
    !>  \compilef{isDarwin}
    !>  \output{isDarwin}
    !>  \include{lineno} example/pm_os/isDarwin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_os](@ref test_pm_os)
    !>
    !>  \final{isDarwin}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isDarwin

    impure module function isDarwin() result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDarwin
#endif
        logical(LK)                     :: itis
    end function

    impure module function isDarwinFailed(failed) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDarwinFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: itis
    end function

    impure module function isDarwinFailedMsg(failed, errmsg) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDarwinFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: itis
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime operating system is Windows.
    !>
    !>  \details
    !>  The identification of an operating system as Windows requires either:<br>
    !>  <ol>
    !>      <li>    the ParaMonte library to have been compiled with the FPP macro `WINDOWS_ENABLED=1`, or,
    !>      <li>    the value contained in the runtime environmental shell variable `OS` match or contain the keyword `windows`.
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `itis`                  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> the runtime operating system is Windows, otherwise, `.false.`.
    !>
    !>  \interface{isWindows}
    !>  \code{.F90}
    !>
    !>      use pm_os, only: isWindows
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: itis
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      itis = isWindows()
    !>      itis = isWindows(failed)
    !>      itis = isWindows(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isLinux](@ref pm_os::isLinux)<br>
    !>  [isDarwin](@ref pm_os::isDarwin)<br>
    !>  [isWindows](@ref pm_os::isWindows)<br>
    !>
    !>  \example{isWindows}
    !>  \include{lineno} example/pm_os/isWindows/main.F90
    !>  \compilef{isWindows}
    !>  \output{isWindows}
    !>  \include{lineno} example/pm_os/isWindows/main.out.F90
    !>
    !>  \test
    !>  [test_pm_os](@ref test_pm_os)
    !>
    !>  \final{isWindows}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isWindows

    impure module function isWindows() result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isWindows
#endif
        logical(LK)                     :: itis
    end function

    impure module function isWindowsFailed(failed) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isWindowsFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: itis
    end function

    impure module function isWindowsFailedMsg(failed, errmsg) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isWindowsFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: itis
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_os ! LCOV_EXCL_LINE