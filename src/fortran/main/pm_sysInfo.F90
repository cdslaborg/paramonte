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
!>  This module contains procedures and generic interfaces for
!>  inferring the operating system kernel type, name, and other information.
!>
!>  \see
!>  [pm_io](@ref pm_io)<br>
!>  [pm_sysInfo](@ref pm_sysInfo)<br>
!>  [pm_sysPath](@ref pm_sysPath)<br>
!>  [pm_sysShell](@ref pm_sysShell)<br>
!>
!>  \test
!>  [test_pm_sysInfo](@ref test_pm_sysInfo)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 3:09 AM, Dec 8, 2017, Dell Medical School, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sysInfo

!>  \cond excluded
!   The Intel FPP defines `linux=1` which conflicts with the `linux` component of `kernelis_type` below.
#if     __INTEL_COMPILER
#undef  linux
#endif
!>  \endcond excluded

    use pm_kind, only: IK, LK, RK, CK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_sysInfo"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [kernelis_type](@ref pm_sysInfo::kernelis_type) class for generating
    !>  objects with `logical` components to determine the operating system (OS) kernel.
    !>
    !>  \details
    !>  When an object of this type is generated, a single object component corresponding to OS kernel is set to `.true.`.<br>
    !>  The OS kernel is inferred either through compile-time preprocessing OS macros if they are defined or otherwise through runtime shell.<br>
    !>
    !>  \warning
    !>  Beware that the environments `"msys"`, `"MinGW"`, `"Cygwin"` are not really standalone Operating Systems,
    !>  rather Linux compatibility environments within Windows platforms.<br>
    !>  As such, when such environments are detected, the `windows`
    !>  component of the object of type [kernelis_type](@ref pm_sysInfo::kernelis_type) is always set to `.true.`.<br>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system kernel.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. It can be present only if `failed` is also present. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `kernelis`              :   The output scalar object of type [kernelis_type](@ref pm_sysInfo::kernelis_type)
    !>                              containing the specifics of the runtime operating system kernel.<br>
    !>
    !>  \interface{kernelis_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysInfo, only: kernelis_type
    !>      character(255, SK) :: errmsg
    !>      type(kernelis_type) :: kernelis
    !>      logical(LK) :: failed
    !>
    !>      kernelis = kernelis_type()
    !>      kernelis = kernelis_type(failed)
    !>      kernelis = kernelis_type(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [getSysInfo](@ref pm_sysInfo::getSysInfo)<br>
    !>  [kernel_type](@ref pm_sysInfo::kernel_type)<br>
    !>  [kernelis_type](@ref pm_sysInfo::kernelis_type)<br>
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [shellis_type](@ref pm_sysShell::shellis_type)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isKernelLinux](@ref pm_sysInfo::isKernelLinux)<br>
    !>  [isKernelDarwin](@ref pm_sysInfo::isKernelDarwin)<br>
    !>  [isKernelWindows](@ref pm_sysInfo::isKernelWindows)<br>
    !>
    !>  \example{kernelis_type}
    !>  \include{lineno} example/pm_sysInfo/kernelis_type/main.F90
    !>  \compilef{kernelis_type}
    !>  \output{kernelis_type}
    !>  \include{lineno} example/pm_sysInfo/kernelis_type/main.out.F90
    !>
    !>  \final{kernelis_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: kernelis_type
        logical(LK) :: windows  = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is Windows.
        logical(LK) :: cygwin   = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is Cygwin.
        logical(LK) :: mingw    = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is MinGW.
        logical(LK) :: msys     = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is MSYS.
        logical(LK) :: linux    = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is Linux.
        logical(LK) :: darwin   = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is Darwin.
        logical(LK) :: freebsd  = .false._LK    !<  \public The scalar `logical` value indicating whether the OS is FreeBSD.
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    interface kernelis_type

    impure module function kernelis_typer() result(kernelis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: kernelis_typer
#endif
        type(kernelis_type)                 :: kernelis
    end function

    impure module function kernelis_typerFailed(failed) result(kernelis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: kernelis_typerFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)       :: failed
        type(kernelis_type)                 :: kernelis
    end function

    impure module function kernelis_typerFailedMsg(failed, errmsg) result(kernelis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: kernelis_typerFailedMsg
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)       :: failed
        character(*,SKG), intent(inout)     :: errmsg
        type(kernelis_type)                 :: kernelis
    end function

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [kernel_type](@ref pm_sysInfo::kernel_type) class for generating
    !>  objects with `logical` components to determine the operating system (OS) kernel and its name.
    !>
    !>  \details
    !>  When an object of this type is generated, a single object component `name` corresponding to OS name is filled with,
    !>  <ol>
    !>      <li>    `"Windows"` if the OS is Windows.
    !>      <li>    `"Cygwin"` if the OS is Cygwin.
    !>      <li>    `"MinGW"` if the OS is MinGW.
    !>      <li>    `"MSYS"` if the OS is MinGW.
    !>      <li>    `"Linux"` if the OS is Linux.
    !>      <li>    `"Darwin"` if the OS is Darwin.
    !>      <li>    `"FreeBSD"` if the OS is FreeBSD.
    !>  </ol>
    !>  Otherwise, if the appropriate preprocessing macro is undefined, the type constructor will
    !>  use the facilities of the runtime shell to infer the Operating System kernel name through either:<br>
    !>  <ol>
    !>      <li>    the Windows `OS` environmental variable in Windows-style runtime shells.
    !>      <li>    the UNIX `OSTYPE` environmental variable in UNIX-style runtime shells.
    !>      <li>    the UNIX `uname` command on unix-style runtime shells.
    !>  </ol>
    !>
    !>  \warning
    !>  The environments `"msys"`, `"MinGW"`, `"Cygwin"` are not really standalone Operating System kernels,
    !>  rather Linux compatibility environments within Windows platforms.<br>
    !>  As such, when such environments are detected, the `windows` component of the `is`
    !>  component of the object of type [kernel_type](@ref pm_sysInfo::kernel_type) is always set to `.true.`.<br>
    !>
    !>  \warning
    !>  Even though Cygwin and MinGW are considered operating systems independently of Windows
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system kernel.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. It can be present only if `failed` is also present. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `kernel`                :   The output scalar object of type [kernel_type](@ref pm_sysInfo::kernel_type)
    !>                              containing the specifics of the operating system kernel name.<br>
    !>
    !>  \interface{kernel_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysInfo, only: kernel_type
    !>      character(255, SK) :: errmsg
    !>      type(kernel_type) :: kernel
    !>      logical(LK) :: failed
    !>
    !>      kernel = kernel_type()
    !>      kernel = kernel_type(failed)
    !>      kernel = kernel_type(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [getSysInfo](@ref pm_sysInfo::getSysInfo)<br>
    !>  [kernel_type](@ref pm_sysInfo::kernel_type)<br>
    !>  [kernelis_type](@ref pm_sysInfo::kernelis_type)<br>
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [shellis_type](@ref pm_sysShell::shellis_type)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isKernelLinux](@ref pm_sysInfo::isKernelLinux)<br>
    !>  [isKernelDarwin](@ref pm_sysInfo::isKernelDarwin)<br>
    !>  [isKernelWindows](@ref pm_sysInfo::isKernelWindows)<br>
    !>
    !>  \example{kernel_type}
    !>  \include{lineno} example/pm_sysInfo/kernel_type/main.F90
    !>  \compilef{kernel_type}
    !>  \output{kernel_type}
    !>  \include{lineno} example/pm_sysInfo/kernel_type/main.out.F90
    !>
    !>  \final{kernel_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type :: kernel_type
        character(:, SK), allocatable       :: name !<  \public The `allocatable` scalar of type `character` of default kind \SK, containing the operating system name.
        type(kernelis_type)                 :: is   !<  \public The scalar object of type [kernelis_type](@ref pm_sysInfo::kernelis_type) containing the logical components indicating the OS.
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    interface kernel_type

    impure module function kernel_typer() result(kernel)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: kernel_typer
#endif
        type(kernel_type)                   :: kernel
    end function

    impure module function kernel_typerFailed(failed) result(kernel)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: kernel_typerFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)       :: failed
        type(kernel_type)                   :: kernel
    end function

    impure module function kernel_typerFailedMsg(failed, errmsg) result(kernel)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: kernel_typerFailedMsg
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)       :: failed
        character(*,SKG), intent(inout)     :: errmsg
        type(kernel_type)                   :: kernel
    end function

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a string containing a comprehensive report of the operating system and platform specifications.
    !>
    !>  \details
    !>  The system information is obtained by first identifying the operating system and
    !>  the runtime shell and then calling one of the following platform/shell-specific commands:<br>
    !>  <ol>
    !>      <li>    The commands `uname -a` + `sysctl -a | grep machdep.cpu` + `system_profiler SPHardwareDataType` on Darwin (macOS) platforms.<br>
    !>      <li>    The commands `uname -a` + `lscpu` + `cat /proc/cpuinfo` on Linux or other Unix-like platforms.
    !>      <li>    The commands `systeminfo` on Windows platforms.
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar of type `logical` of default kind \LK, that is `.true.` if and only if the process of
    !>                              fetching system information fails.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. It can be present **only if** the input argument `failed` is also present.)
    !>
    !>  \return
    !>  `sysInfo`               :   The output `allocatable` scalar `character` of default kind \SK containing the system information.<br>
    !>
    !>  \interface{getSysInfo}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use test_pm_sys, only: getSysInfo
    !>
    !>      character(:, SK), allocatable :: sysInfo
    !>
    !>      sysInfo = getSysInfo()
    !>      sysInfo = getSysInfo(failed)
    !>      sysInfo = getSysInfo(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getSysInfo](@ref pm_sysInfo::getSysInfo)<br>
    !>  [isKernelLinux](@ref pm_sysInfo::isKernelLinux)<br>
    !>  [isKernelDarwin](@ref pm_sysInfo::isKernelDarwin)<br>
    !>  [isKernelWindows](@ref pm_sysInfo::isKernelWindows)<br>
    !>  [kernelis_type](@ref pm_sysInfo::kernelis_type)<br>
    !>  [kernel_type](@ref pm_sysInfo::kernel_type)<br>
    !>
    !>  \example{getSysInfo}
    !>  \include{lineno} example/pm_sysInfo/getSysInfo/main.F90
    !>  \compilef{getSysInfo}
    !>  \output{getSysInfo}
    !>  \include{lineno} example/pm_sysInfo/getSysInfo/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sys](@ref test_pm_sys)
    !>
    !>  \final{getSysInfo}
    !>
    !>  \author
    !>  Amir Shahmoradi, Tuesday March 7, 2017, 3:09 AM, ICES, The University of Texas at Austin

    interface getSysInfo

    module function getSysInfo() result(sysInfo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSysInfo
#endif
        use pm_kind, only: SKG => SK
        character(:,SKG)                , allocatable   :: sysInfo
    end function

    module function getSysInfoFailed(failed) result(sysInfo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSysInfoFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)                   :: failed
        character(:,SKG)                , allocatable   :: sysInfo
    end function

    module function getSysInfoFailedMsg(failed, errmsg) result(sysInfo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSysInfoFailedMsg
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)                   :: failed
        character(*,SKG), intent(inout)                 :: errmsg
        character(:,SKG)                , allocatable   :: sysInfo
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the Operating System kernel is the Windows.
    !>
    !>  \details
    !>  Refer to the documentation details of [kernel_type](@ref pm_sysInfo::kernel_type)
    !>  for information on the operating system kernels that recognized and supported by the ParaMonte library.<br>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system kernel type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `itis`                  :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the operating system kernel is Windows-based.
    !>
    !>  \interface{isKernelWindows}
    !>  \code{.F90}
    !>
    !>      use pm_sysInfo, only: isKernelWindows
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: itis
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      itis = isKernelWindows()
    !>      itis = isKernelWindows(failed)
    !>      itis = isKernelWindows(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getSysInfo](@ref pm_sysInfo::getSysInfo)<br>
    !>  [isKernelLinux](@ref pm_sysInfo::isKernelLinux)<br>
    !>  [isKernelDarwin](@ref pm_sysInfo::isKernelDarwin)<br>
    !>  [isKernelWindows](@ref pm_sysInfo::isKernelWindows)<br>
    !>  [kernelis_type](@ref pm_sysInfo::kernelis_type)<br>
    !>  [kernel_type](@ref pm_sysInfo::kernel_type)<br>
    !>
    !>  \example{isKernelWindows}
    !>  \include{lineno} example/pm_sysInfo/isKernelWindows/main.F90
    !>  \compilef{isKernelWindows}
    !>  \output{isKernelWindows}
    !>  \include{lineno} example/pm_sysInfo/isKernelWindows/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysInfo](@ref test_pm_sysInfo)
    !>
    !>  \final{isKernelWindows}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface isKernelWindows

    module function isKernelWindows() result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelWindows
#endif
        logical(LK)                     :: itis
    end function

    module function isKernelWindowsFailed(failed) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelWindowsFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: itis
    end function

    module function isKernelWindowsFailedMsg(failed, errmsg) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelWindowsFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: itis
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the Operating System kernel is the macOS Darwin.
    !>
    !>  \details
    !>  Refer to the documentation details of [kernel_type](@ref pm_sysInfo::kernel_type)
    !>  for information on the operating system kernels that recognized and supported by the ParaMonte library.<br>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system kernel type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `itis`                  :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the operating system kernel is Darwin-based.
    !>
    !>  \interface{isKernelDarwin}
    !>  \code{.F90}
    !>
    !>      use pm_sysInfo, only: isKernelDarwin
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: itis
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      itis = isKernelDarwin()
    !>      itis = isKernelDarwin(failed)
    !>      itis = isKernelDarwin(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getSysInfo](@ref pm_sysInfo::getSysInfo)<br>
    !>  [isKernelLinux](@ref pm_sysInfo::isKernelLinux)<br>
    !>  [isKernelDarwin](@ref pm_sysInfo::isKernelDarwin)<br>
    !>  [isKernelWindows](@ref pm_sysInfo::isKernelWindows)<br>
    !>  [kernelis_type](@ref pm_sysInfo::kernelis_type)<br>
    !>  [kernel_type](@ref pm_sysInfo::kernel_type)<br>
    !>
    !>  \example{isKernelDarwin}
    !>  \include{lineno} example/pm_sysInfo/isKernelDarwin/main.F90
    !>  \compilef{isKernelDarwin}
    !>  \output{isKernelDarwin}
    !>  \include{lineno} example/pm_sysInfo/isKernelDarwin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysInfo](@ref test_pm_sysInfo)
    !>
    !>  \final{isKernelDarwin}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface isKernelDarwin

    module function isKernelDarwin() result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelDarwin
#endif
        logical(LK)                     :: itis
    end function

    module function isKernelDarwinFailed(failed) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelDarwinFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: itis
    end function

    module function isKernelDarwinFailedMsg(failed, errmsg) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelDarwinFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: itis
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the Operating System kernel is Linux.
    !>
    !>  \details
    !>  Refer to the documentation details of [kernel_type](@ref pm_sysInfo::kernel_type)
    !>  for information on the operating system kernels that recognized and supported by the ParaMonte library.<br>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the operating system kernel type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `itis`         :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the operating system kernel is Linux-based.
    !>
    !>  \interface{isKernelLinux}
    !>  \code{.F90}
    !>
    !>      use pm_sysInfo, only: isKernelLinux
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: itis
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      itis = isKernelLinux()
    !>      itis = isKernelLinux(failed)
    !>      itis = isKernelLinux(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getSysInfo](@ref pm_sysInfo::getSysInfo)<br>
    !>  [isKernelLinux](@ref pm_sysInfo::isKernelLinux)<br>
    !>  [isKernelDarwin](@ref pm_sysInfo::isKernelDarwin)<br>
    !>  [isKernelWindows](@ref pm_sysInfo::isKernelWindows)<br>
    !>  [kernelis_type](@ref pm_sysInfo::kernelis_type)<br>
    !>  [kernel_type](@ref pm_sysInfo::kernel_type)<br>
    !>
    !>  \example{isKernelLinux}
    !>  \include{lineno} example/pm_sysInfo/isKernelLinux/main.F90
    !>  \compilef{isKernelLinux}
    !>  \output{isKernelLinux}
    !>  \include{lineno} example/pm_sysInfo/isKernelLinux/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysInfo](@ref test_pm_sysInfo)
    !>
    !>  \final{isKernelLinux}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    interface isKernelLinux

    module function isKernelLinux() result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelLinux
#endif
        logical(LK)                     :: itis
    end function

    module function isKernelLinuxFailed(failed) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelLinuxFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: itis
    end function

    module function isKernelLinuxFailedMsg(failed, errmsg) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isKernelLinuxFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: itis
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
#undef PURITY
!>  \endcond excluded

end module pm_sysInfo ! LCOV_EXCL_LINE