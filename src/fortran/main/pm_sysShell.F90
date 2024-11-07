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
!>  inferring the runtime system shell type and fetching information from the shell.
!>
!>  \see
!>  [pm_os](@ref pm_os)<br>
!>  [pm_io](@ref pm_io)<br>
!>  [pm_sysInfo](@ref pm_sysInfo)<br>
!>  [pm_sysPath](@ref pm_sysPath)<br>
!>  [pm_sysShell](@ref pm_sysShell)<br>
!>
!>  \test
!>  [test_pm_sysShell](@ref test_pm_sysShell)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 3:09 AM, Dec 8, 2017, Dell Medical School, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sysShell

    use pm_kind, only: IK, LK, RK, CK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_sysShell"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant vector of scalars of type `character` of default kind \SK,
    !>  containing the Windows shell environmental variable names that potentially
    !>  the path to the system directory within which temporary files can be stored.
    !>
    !>  \devnote
    !>  The specific order of the vector elements is important and signifies the popularity of the environment variable names.<br>
    !>
    !>  \see
    !>  [VARENV_DIRTEMP_UNIX](@ref pm_sysShell::VARENV_DIRTEMP_UNIX)<br>
    !>  [VARENV_DIRTEMP_WINDOWS](@ref pm_sysShell::VARENV_DIRTEMP_WINDOWS)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedGetEnvVar](@ref pm_sysShell::isFailedGetEnvVar)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isFailedExec](@ref pm_sysShell::isFailedExec)<br>
    !>  [getPathTemp](@ref pm_sysPath::getPathTemp)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>
    !>  \final{VARENV_DIRTEMP_WINDOWS}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(4, SK), parameter :: VARENV_DIRTEMP_WINDOWS(*) = [character(4, SK) :: "TEMP", "TMP"]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: VARENV_DIRTEMP_WINDOWS
#endif

    !>  \brief
    !>  The constant vector of scalars of type `character` of default kind \SK,
    !>  containing the Unix shell environmental variable names that potentially
    !>  the path to the system directory within which temporary files can be stored.
    !>  \devnote
    !>  The specific order of the vector elements is important and signifies the popularity of the environment variable names.<br>
    !>
    !>  \see
    !>  [VARENV_DIRTEMP_UNIX](@ref pm_sysShell::VARENV_DIRTEMP_UNIX)<br>
    !>  [VARENV_DIRTEMP_WINDOWS](@ref pm_sysShell::VARENV_DIRTEMP_WINDOWS)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedGetEnvVar](@ref pm_sysShell::isFailedGetEnvVar)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isFailedExec](@ref pm_sysShell::isFailedExec)<br>
    !>  [getPathTemp](@ref pm_sysPath::getPathTemp)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>
    !>  \final{VARENV_DIRTEMP_UNIX}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(7, SK), parameter :: VARENV_DIRTEMP_UNIX(*) = [character(7, SK) :: "TMPDIR", "TMP", "TEMP", "TEMPDIR"]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: VARENV_DIRTEMP_UNIX
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [shellis_type](@ref pm_sysShell::shellis_type) class for generating
    !>  objects with `logical` components to determine the runtime shell type of the operating system.
    !>
    !>  \details
    !>  When an object of this type is generated, a single object component corresponding to runtime shell name is set to `.true.`.<br>
    !>  If the runtime shell name is determined to be `sh`, additional attempt will be made to resolve the actual runtime shell name.<br>
    !>  This is because `sh` is frequently a symlink with an actual target (for example `dash` on Ubuntu).<br>
    !>
    !>  POSIX-compliant shells
    !>  ----------------------
    !>
    !>  Additionally, the `posix` component of the object is set to `.true.` if the runtime shell name corresponds to,
    !>  <ol>
    !>      <li>    ash
    !>      <li>    bash
    !>      <li>    csh
    !>      <li>    dash
    !>      <li>    ksh
    !>      <li>    pwsh (Microsoft PowerShell core)
    !>      <li>    sh
    !>      <li>    tcsh
    !>      <li>    zsh
    !>      <li>    yash
    !>  </ol>
    !>  Normally (, but not necessarily), the runtime operating system for the above shells is Unix (either macOS or variants of linux).<br>
    !>  As of July 2022, the above shells are known to have the highest level of compliance with the
    !>  [POSIX](https://pubs.opengroup.org/onlinepubs/9699919799/utilities/V3_chap01.html) standard.<br>
    !>
    !>  Windows-compliant shells
    !>  ------------------------
    !>
    !>  Alternatively, the `windows` component of the object is set to `.true.` if the runtime shell corresponds to,
    !>  <ol>
    !>      <li>    CMD
    !>      <li>    Windows PowerShell
    !>  </ol>
    !>  As of July 2022, the above shells are exclusively available on Windows operating systems.<br>
    !>
    !>  Other shells
    !>  ------------
    !>
    !>  As of July 2022, shells that adhere to neither Windows nor POSIX conventions are the following,
    !>  <ol>
    !>      <li>    fish
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell name.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. It can be present only if `failed` is also present. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `shell`                 :   The output scalar object of type [shellis_type](@ref pm_sysShell::shellis_type) containing the specifics of the runtime system shell.<br>
    !>
    !>  \interface{shellis_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysShell, only: shellis_type
    !>      character(255, SK) :: errmsg
    !>      type(shellis_type) :: shell
    !>      logical(LK) :: failed
    !>
    !>      shell = shellis_type()
    !>      shell = shellis_type(failed)
    !>      shell = shellis_type(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [Comparison of command shells](https://en.wikipedia.org/wiki/Comparison_of_command_shells)<br>
    !>
    !>  \example{shellis_type}
    !>  \include{lineno} example/pm_sysShell/shellis_type/main.F90
    !>  \compilef{shellis_type}
    !>  \output{shellis_type}
    !>  \include{lineno} example/pm_sysShell/shellis_type/main.out.F90
    !>
    !>  \final{shellis_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: shellis_type
        logical(LK)                     :: ash        = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Almquist shell (**dash**).
        logical(LK)                     :: bash       = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix GNU Bourne-Again shell (**Bash**).
        logical(LK)                     :: cmd        = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Windows Command Prompt (**CMD**).
        logical(LK)                     :: csh        = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix C shell (**csh**).
        logical(LK)                     :: dash       = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Debian Almquist shell (**dash**).
        logical(LK)                     :: fish       = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Friendly Interactive SHell (**fish**).
        logical(LK)                     :: ksh        = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Korn shell (**ksh**).
        logical(LK)                     :: posix      = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is POSIX-compliant (that is, Unix-like, partially or fully POSIX-compliant, e.g., pwsh (PowerShell Core), sh, bash, csh, zsh, fish, ...).
        logical(LK)                     :: powershell = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is the Microsoft PowerShell (**pwsh**).
        logical(LK)                     :: sh         = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Bourne shell (**sh**).
                                                                      !!          \note
                                                                      !!          The Shell Command Language (`sh`) is a programming language described by the POSIX standard.<br>
                                                                      !!          Because `sh` is a specification, not an implementation, `/bin/sh` is a symlink (or a hard link) to an actual implementation on most POSIX systems.<br>
                                                                      !!          The value `sh = .true.` implies that the runtime shell most likely points to a POSIX-compliant shell. To find out the target of `sh`, use the command `file -h /bin/sh`.<br>
        logical(LK)                     :: tcsh       = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix `tcsh` shell (**tcsh**).
        logical(LK)                     :: windows    = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Windows-compliant (e.g., pwsh on Windows (Windows PowerShell), CMD).
        logical(LK)                     :: zsh        = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Z shell (**zsh**).
        logical(LK)                     :: yash       = .false._LK    !<  \public The scalar `logical` value indicating whether the shell is Unix Yet Another Shell (**yash**).
    end type

    !>  \cond excluded
    interface shellis_type

    impure module function shellis_typer() result(shellis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shellis_typer
#endif
        type(shellis_type)              :: shellis
    end function

    impure module function shellis_typerFailed(failed) result(shellis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shellis_typerFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        type(shellis_type)              :: shellis
    end function

    impure module function shellis_typerFailedMsg(failed, errmsg) result(shellis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shellis_typerFailedMsg
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        character(*,SKG), intent(inout) :: errmsg
        type(shellis_type)              :: shellis
    end function

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [shell_type](@ref pm_sysShell::shell_type) class for
    !>  generating objects to determine the runtime shell type of the operating system.
    !>
    !>  \details
    !>  Note that the system shells are frequently cross-platform.<br>
    !>  For example,
    !>  <ol>
    !>      <li>    Bash shell is frequently also available on Windows platforms.
    !>      <li>    Microsoft PowerShell Core (pwsh) is compliant with **IEEE POSIX 1003.2 standard for Unix shells** and can be installed on all major operating systems.
    !>  </ol>
    !>
    !>  See the documentation of [shellis_type](@ref pm_sysShell::shellis_type) for details and meaning of the `is` component of objects of this type.<br>
    !>  Note that both gfortran as of version 12 and Intel ifort as of version 2022 use the Windows Command prompt as the default shell,
    !>  even if the program is compiled and run from within a POSIX-compliant shell (like Git Bash).
    !>
    !>  Frequently, the runtime shell `sh` is merely a symlink to a more modern shell compatible with Bourne shell.<br>
    !>  To ensure the target shell name is correctly inferred, the program will report the target of `sh` if it is a symlink.<br>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell name.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. It can be present only if `failed` is also present. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `shell`                 :   The output scalar object of type [shell_type](@ref pm_sysShell::shell_type) containing the specifics of the runtime system shell.<br>
    !>
    !>  \interface{shell_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysShell, only: shell_type
    !>      character(255, SK) :: errmsg
    !>      type(shell_type) :: shell
    !>      logical(LK) :: failed
    !>
    !>      shell = shell_type()
    !>      shell = shell_type(failed)
    !>      shell = shell_type(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  Inferring the shell type at runtime is computationally expensive. Therefore, to avoid this extra computational cost,
    !>  the procedures under this generic interface construct a static object of type [shell_type](@ref pm_sysShell::shell_type)
    !>  which is `save`d within the procedures for internal usage in subsequent calls to the procedures.<br>
    !>  The presence of such a sticky internal object should be fine even in threaded (e.g., OpenMP) applications,
    !>  since the state of the object is initialized once and not changed anymore throughout the entire program life.
    !>
    !>  \see
    !>  [shellis_type](@ref pm_sysShell::shellis_type)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [Comparison of command shells](https://en.wikipedia.org/wiki/Comparison_of_command_shells)<br>
    !>
    !>  \example{shell_type}
    !>  \include{lineno} example/pm_sysShell/shell_type/main.F90
    !>  \compilef{shell_type}
    !>  \output{shell_type}
    !>  \include{lineno} example/pm_sysShell/shell_type/main.out.F90
    !>
    !>  \final{shell_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: shell_type
        type(shellis_type)              :: is               !<  \public The scalar object of type [shellis_type](@ref pm_sysShell::shellis_type) whose `logical` components indicate the shell name.
        character(1, SK)                :: dirsep = SK_"/"  !<  \public The scalar `character` of default kind \SK of length `1` containing the **preferred** single-character shell **directory separator** (e.g., `\` (Windows) or `/` (POSIX or if the shell type is unknown)).
        character(1, SK)                :: pathsep = SK_":" !<  \public The scalar `character` of default kind \SK of length `1` containing the **preferred** single-character shell **path separator** (e.g., `;` (Windows) or `:` (POSIX or if the shell type is unknown)).
        character(:, SK), allocatable   :: dirseps          !<  \public The `allocatable` scalar `character` of default kind \SK of arbitrary length containing all supported single-character shell **directory separators** (e.g., `\/` (Windows) or `/` (POSIX or if the shell type is unknown)).
        character(:, SK), allocatable   :: name             !<  \public The `allocatable` scalar `character` containing the name of or path to the current shell.
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    interface shell_type

    impure module function shell_typer() result(shell)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shell_typer
#endif
        type(shell_type)                :: shell
    end function

    impure module function shell_typerFailed(failed) result(shell)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shell_typerFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        type(shell_type)                :: shell
    end function

    impure module function shell_typerFailedMsg(failed, errmsg) result(shell)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shell_typerFailedMsg
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        character(*,SKG), intent(inout) :: errmsg
        type(shell_type)                :: shell
    end function

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the value of the input environment variable `name`
    !>  can be successfully retrieved as a trimmed string in the output argument `value`.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.<br>
    !>
    !>  \details
    !>  The functionality of the procedures under this generic interface are highly similar to that of the Fortran intrinsic `get_environment_variable()`.<br>
    !>  However, the procedures of this generic interface remove the burden of guessing the correct size of the output string value by try and error.<br>
    !>
    !>  \param[in]      name    :   The input scalar `character` of default kind \SK containing the name of the
    !>                              environment variable whose value should be retrieved from the runtime system shell.<br>
    !>                              The argument `name` must neither be empty nor all blank characters.
    !>  \param[out]     value   :   The output `allocatable` scalar `character` of default kind \SK containing the value of the environment variable `name`.<br>
    !>                              The output `value` is guaranteed to capture the full value of the environment variable `name` and without any trailing blanks.<br>
    !>                              **The `allocation` status or the contents of `value` is undefined if a runtime error occurs.**<br>
    !>                              **Always check the returned value of `failed` before using `value`.**
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. It can be present **only if** `failed` is also present.)
    !>  \param[in]      length  :   The input scalar `integer` of default kind \IK representing the best-guess length of the output `value` of the environment variable.<br>
    !>                              By default, the initial length of `value` is `2**13 - 1 = 8197` characters. This is suitable for very long environment variables like `PATH`.<br>
    !>                              However, this may be too large in many other cases. If so, the user can provide a better starting guess for the length of `value`
    !>                              by setting the input argument `length` to a proper best guess positive `integer`.<br>
    !>                              (**optional**, default = `8191`)
    !>
    !>  \return
    !>  `failed`                :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> the process of retrieving the value of the environment variable fails at any point.
    !>
    !>  \interface{isFailedGetEnvVar}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedGetEnvVar
    !>      character(:, SK), allocatable :: value
    !>      character(255, SK) :: errmsg
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedGetEnvVar(name, value, length = length)
    !>      failed = isFailedGetEnvVar(name, value, errmsg, length = length)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `length > 0_IK` must hold.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isFailedPutEnvVar](@ref pm_sysShell::isFailedPutEnvVar)<br>
    !>
    !>  \example{isFailedGetEnvVar}
    !>  \include{lineno} example/pm_sysShell/isFailedGetEnvVar/main.F90
    !>  \compilef{isFailedGetEnvVar}
    !>  \output{isFailedGetEnvVar}
    !>  \include{lineno} example/pm_sysShell/isFailedGetEnvVar/main.out.F90
    !>
    !>  \final{isFailedGetEnvVar}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGetEnvVar

    impure module function isFailedGetEnvVar(name, value, length) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetEnvVar
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                    :: name
        character(:,SKG), intent(out)   , allocatable   :: value
        integer(IK)     , intent(in)    , optional      :: length
        logical(LK)                                     :: failed
    end function

    impure module function isFailedGetEnvVarMsg(name, value, errmsg, length) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetEnvVarMsg
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                    :: name
        character(:,SKG), intent(out)   , allocatable   :: value
        character(*, SK), intent(inout)                 :: errmsg
        integer(IK)     , intent(in)    , optional      :: length
        logical(LK)                                     :: failed
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the attempt to (re)define the environment variable `name`
    !>  with the specified `value` in the runtime shell **fails, otherwise return `.false.` upon success.<br>
    !>
    !>  \param[in]      name    :   The input scalar `character` of default kind \SK containing the name of the
    !>                              environment variable whose value should be retrieved from the runtime system shell.<br>
    !>                              The argument `name` must neither be empty nor all blank characters.
    !>  \param[out]     value   :   The output `allocatable` scalar `character` of default kind \SK containing the value of the environment variable `name`.<br>
    !>                              The output `value` is guaranteed to capture the full value of the environment variable `name` and without any trailing blanks.<br>
    !>                              **If the environment variable does not exist or a runtime error occurs, `value` will be allocated to empty string upon return.**
    !>
    !>  \return
    !>  `failed`                :   The output scalar `logical` of default kind \LK.<br>
    !>                              It is `.true.` <b>if and only if</b> an error occurs while defining
    !>                              the environment variable `name` with the specified `value` in the runtime shell.<br>
    !>
    !>  \interface{isFailedPutEnvVar}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedPutEnvVar
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedPutEnvVar(name, value)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isFailedGetEnvVar](@ref pm_sysShell::isFailedGetEnvVar)<br>
    !>
    !>  \example{isFailedPutEnvVar}
    !>  \include{lineno} example/pm_sysShell/isFailedPutEnvVar/main.F90
    !>  \compilef{isFailedPutEnvVar}
    !>  \output{isFailedPutEnvVar}
    !>  \include{lineno} example/pm_sysShell/isFailedPutEnvVar/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \todo
    !>  \plow
    !>  A subroutine version of this functional interface could be implemented in
    !>  future to avoid allocations and allow for non-default character kinds.<br>
    !>
    !>  \final{isFailedPutEnvVar}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedPutEnvVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isFailedPutEnvVar(name, value) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedPutEnvVar
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: name, value
        logical(LK)                     :: failed
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the output of the input system command is successfully retrieved in the `allocatable` argument `output`.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.
    !>
    !>  \details
    !>  To do so, the procedures under this generic interface call the system to run the requested command
    !>  while redirecting its output to a temporary file for subsequent retrieval and return to the user.
    !>
    !>  \warning
    !>  **Use this functionality responsibly and with caution.**<br>
    !>  Due to their simplicity, the procedures under this generic interface
    !>  have high potential of abuse by passing system commands that can wreak havoc.<br>
    !>  The responsible usage of this functionality must be currently guaranteed by the user.
    !>
    !>  \param[in]      command     :   The input scalar `character` of default kind \SK containing the system command whose output must be retrieved.
    !>  \param[out]     output      :   The output `allocatable` scalar `character` of default kind \SK containing the output from the system `command` call.<br>
    !>                                  The `output` may contain the text redirected to the standard error if the specified `command` fails, even though the processor calls it and returns successfully.<br>
    !>                                  This can be verified by testing whether the `optional` output argument `exitstat` is non-zero (implying the failure of the specified `command` to return successfully).<br>
    !>                                  The allocation status of `output` remains undefined if the procedure returns `failed` set to `.true.`.<br>
    !>                                  **Always check the returned value of `failed` before using `output`.**<br>
    !>                                  This behavior is sensible because when the procedure fails, the contents of `output` is garbage and should not be subsequently used.<br>
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                                  A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output if an error occurs.)
    !>  \param[out]     exitstat    :   The output scalar `integer` of default kind \IK, representing the exit code returned by the input `command`.<br>
    !>                                  For more information, see the corresponding argument of [isFailedExec](@ref pm_sysShell::isFailedExec).<br>
    !>                                  The specific returned value of `exitstat` is processor and application dependent.<br>
    !>                                  However, **a zero `exitstat` frequently implies the specified `command` returned the control successfully to the processor, while non-zero `exitstat` implies `command` failure**.<br>
    !>                                  See also the corresponding output argument of [isFailedExec](@ref pm_sysShell::isFailedExec).<br>
    !>                                  (**optional**, If missing, user may find it difficult to tell whether the contents of `output` are from stdout or strerr.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> the process of retrieving the command output fails at any point.<br>
    !>
    !>  \interface{isFailedGetOutput}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedGetOutput
    !>
    !>      failed = isFailedGetOutput(command, output, exitstat = exitstat)
    !>      failed = isFailedGetOutput(command, output, errmsg, exitstat = exitstat)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  User must always verify the function output `failed` is `.false.` before using the `output`.
    !>
    !>  \example{isFailedGetOutput}
    !>  \include{lineno} example/pm_sysShell/isFailedGetOutput/main.F90
    !>  \compilef{isFailedGetOutput}
    !>  \output{isFailedGetOutput}
    !>  \include{lineno} example/pm_sysShell/isFailedGetOutput/main.out.F90
    !>
    !>  \final{isFailedGetOutput}
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface should be extended to support character `output` of non-default kind once non-default characters is support by all relevant Fortran compilers.
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGetOutput

    impure module function isFailedGetOutput(command, output, exitstat) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetOutput
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                    :: command
        character(:,SKG), intent(out)   , allocatable   :: output
        integer(IK)     , intent(out)   , optional      :: exitstat
        logical(LK)                                     :: failed
    end function

    impure module function isFailedGetOutputMsg(command, output, errmsg, exitstat) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetOutputMsg
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                    :: command
        character(:,SKG), intent(out)   , allocatable   :: output
        character(*, SK), intent(inout)                 :: errmsg
        integer(IK)     , intent(out)   , optional      :: exitstat
        logical(LK)                                     :: failed
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the path to the **shell-specified temporary directory** is successfully returned in the output `dirTemp` argument.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.
    !>
    !>  \details
    !>  The **shell-specified temporary directory** returned by the procedures under this generic interface is suitable for storing the temporary files and folders.<br>
    !>  The procedures under this generic interface return the path **designated** by the current runtime system shell as temporary folder.<br>
    !>  Therefore, the temporary directory returned by the procedures under this generic interface is platform-dependent.<br>
    !>  -#  In POSIX shells, it is the directory specified in one of the following the environment variables `TMPDIR`, `TMP`, `TEMP`, `TEMPDIR`.<br>
    !>      If none of the environment variables exist, the default temporary directory `"/tmp"` is returned as the temporary directory **only if it exists**.<br>
    !>      Otherwise, the procedure returns `failed = .true.` and `dirTemp` variable remains unallocated.<br>
    !>  -#  In Windows shells, it is the directory stored in the environment variable `TEMP` or `TMP`.<br>
    !>      If none of the environment variables exist, the procedure returns `failed = .true.` and `dirTemp` variable remains unallocated.
    !>
    !>  \param[out]     dirTemp :   The output `allocatable` scalar `character` of default kind \SK containing
    !>                              the path to the temporary directory specified by the runtime system shell.<br>
    !>                              **The `allocation` status or the value of `dirTemp` is undefined if a runtime error occurs.**<br>
    !>                              **Always check the returned value of `failed` before using `dirTemp`.**
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `failed`                :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> the process of retrieving the system shell temporary path fails at any point.<br>
    !>
    !>  \interface{isFailedGetDirTemp}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedGetDirTemp
    !>      use pm_kind, only: LK
    !>      character(:, SK), allocatable :: dirTemp
    !>      character(255, SK) :: errmsg
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedGetDirTemp(dirTemp)
    !>      failed = isFailedGetDirTemp(dirTemp, errmsg)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  User must always verify the function output `failed` is `.false.` before using the output `dirTemp`.
    !>
    !>  \impure
    !>
    !>  \see
    !>  [VARENV_DIRTEMP_UNIX](@ref pm_sysShell::VARENV_DIRTEMP_UNIX)<br>
    !>  [VARENV_DIRTEMP_WINDOWS](@ref pm_sysShell::VARENV_DIRTEMP_WINDOWS)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedGetEnvVar](@ref pm_sysShell::isFailedGetEnvVar)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>  [isShellPosix](@ref pm_sysShell::isShellPosix)<br>
    !>  [isFailedExec](@ref pm_sysShell::isFailedExec)<br>
    !>  [getPathTemp](@ref pm_sysPath::getPathTemp)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>
    !>  \example{isFailedGetDirTemp}
    !>  \include{lineno} example/pm_sysShell/isFailedGetDirTemp/main.F90
    !>  \compilef{isFailedGetDirTemp}
    !>  \output{isFailedGetDirTemp}
    !>  \include{lineno} example/pm_sysShell/isFailedGetDirTemp/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \todo
    !>  \plow
    !>  A subroutine version of this functional interface could be implemented in the future to avoid allocations and
    !>  allow for non-default character kinds (once the Fortran compilers fully support non-default characters).
    !>
    !>  \final{isFailedGetDirTemp}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGetDirTemp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isFailedGetDirTemp(dirTemp) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetDirTemp
#endif
        use pm_kind, only: SKG => SK
        character(:,SKG), intent(out)   , allocatable   :: dirTemp
        logical(LK)                                     :: failed
    end function

    impure module function isFailedGetDirTempMsg(dirTemp, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetDirTempMsg
#endif
        use pm_kind, only: SKG => SK
        character(:,SKG), intent(out)   , allocatable   :: dirTemp
        character(*, SK), intent(inout)                 :: errmsg
        logical(LK)                                     :: failed
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime system shell is POSIX-compliant **regardless** of the Operating System.
    !>
    !>  \details
    !>  Refer to the documentation details of [shell_type](@ref pm_sysShell::shell_type)
    !>  for information on what system shells are considered POSIX-compliant on what operating systems.<br>
    !>  The procedures of this generic interface serve as a shortcut to the `posix` component of [shellis_type](@ref pm_sysShell::shellis_type).
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `shellIsPosix`          :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the runtime system shell is POSIX-based (e.g., PowerShell Core, sh, bash, chs, dash, zsh, ...).
    !>
    !>  \interface{isShellPosix}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isShellPosix
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: shellIsPosix
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      shellIsPosix = isShellPosix()
    !>      shellIsPosix = isShellPosix(failed)
    !>      shellIsPosix = isShellPosix(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>
    !>  \example{isShellPosix}
    !>  \include{lineno} example/pm_sysShell/isShellPosix/main.F90
    !>  \compilef{isShellPosix}
    !>  \output{isShellPosix}
    !>  \include{lineno} example/pm_sysShell/isShellPosix/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \final{isShellPosix}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isShellPosix

    impure module function isShellPosix() result(shellIsPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellPosix
#endif
        logical(LK)                     :: shellIsPosix
    end function

    impure module function isShellPosixFailed(failed) result(shellIsPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellPosixFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: shellIsPosix
    end function

    impure module function isShellPosixFailedMsg(failed, errmsg) result(shellIsPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellPosixFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: shellIsPosix
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime system shell is Windows-based (e.g., CMD, PowerShell) regardless of the Operating System.
    !>
    !>  \details
    !>  Refer to the documentation details of [shell_type](@ref pm_sysShell::shell_type)
    !>  for information on what system shells are considered Windows-compliant on what operating systems.<br>
    !>  The procedures of this generic interface serve as a shortcut to the `windows` component of [shellis_type](@ref pm_sysShell::shellis_type).
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `shellIsWindows`        :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the runtime system shell is Windows-based (e.g., CMD, Windows PowerShell).
    !>
    !>  \interface{isShellWindows}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isShellWindows
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: shellIsWindows
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      shellIsWindows = isShellWindows()
    !>      shellIsWindows = isShellWindows(failed)
    !>      shellIsWindows = isShellWindows(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>
    !>  \example{isShellWindows}
    !>  \include{lineno} example/pm_sysShell/isShellWindows/main.F90
    !>  \compilef{isShellWindows}
    !>  \output{isShellWindows}
    !>  \include{lineno} example/pm_sysShell/isShellWindows/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \final{isShellWindows}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isShellWindows

    impure module function isShellWindows() result(shellIsWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellWindows
#endif
        logical(LK)                     :: shellIsWindows
    end function

    impure module function isShellWindowsFailed(failed) result(shellIsWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellWindowsFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: shellIsWindows
    end function

    impure module function isShellWindowsFailedMsg(failed, errmsg) result(shellIsWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellWindowsFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: shellIsWindows
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime system shell is Windows CMD.
    !>
    !>  \details
    !>  Refer to the documentation details of [shell_type](@ref pm_sysShell::shell_type)
    !>  for information on what system shells are available on what operating systems.<br>
    !>  The procedures of this generic interface serve as a shortcut to the `cmd` component of [shellis_type](@ref pm_sysShell::shellis_type).
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `shellIsCMD`            :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the runtime system shell is Windows-based (e.g., CMD, PowerShell).
    !>
    !>  \interface{isShellCMD}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isShellCMD
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: shellIsCMD
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      shellIsCMD = isShellCMD()
    !>      shellIsCMD = isShellCMD(failed)
    !>      shellIsCMD = isShellCMD(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>
    !>  \example{isShellCMD}
    !>  \include{lineno} example/pm_sysShell/isShellCMD/main.F90
    !>  \compilef{isShellCMD}
    !>  \output{isShellCMD}
    !>  \include{lineno} example/pm_sysShell/isShellCMD/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \final{isShellCMD}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isShellCMD

    impure module function isShellCMD() result(shellIsCMD)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellCMD
#endif
        logical(LK)                     :: shellIsCMD
    end function

    impure module function isShellCMDFailed(failed) result(shellIsCMD)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellCMDFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: shellIsCMD
    end function

    impure module function isShellCMDFailedMsg(failed, errmsg) result(shellIsCMD)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellCMDFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: shellIsCMD
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the runtime system shell is Microsoft PowerShell Core (on Unix systems) or Windows PowerShell (on Windows systems).
    !>
    !>  \details
    !>  Refer to the documentation details of [shell_type](@ref pm_sysShell::shell_type)
    !>  for information on what system shells are available on what operating systems.<br>
    !>  The procedures of this generic interface serve as a shortcut to the `powershell` component of [shellis_type](@ref pm_sysShell::shellis_type).
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell type.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `shellIsPowerShell`     :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b>
    !>                              the runtime system shell is either PowerShell Core (on Unix) or Windows PowerShell (on Windows).
    !>
    !>  \interface{isShellPowerShell}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isShellPowerShell
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>      logical(LK) :: shellIsPowerShell
    !>      character(255, SK) :: errmsg = SK_""
    !>
    !>      shellIsPowerShell = isShellPowerShell()
    !>      shellIsPowerShell = isShellPowerShell(failed)
    !>      shellIsPowerShell = isShellPowerShell(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [isShellWindows](@ref pm_sysShell::isShellWindows)<br>
    !>
    !>  \example{isShellPowerShell}
    !>  \include{lineno} example/pm_sysShell/isShellPowerShell/main.F90
    !>  \compilef{isShellPowerShell}
    !>  \output{isShellPowerShell}
    !>  \include{lineno} example/pm_sysShell/isShellPowerShell/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \final{isShellPowerShell}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isShellPowerShell

    impure module function isShellPowerShell() result(shellIsPowerShell)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellPowerShell
#endif
        logical(LK)                     :: shellIsPowerShell
    end function

    impure module function isShellPowerShellFailed(failed) result(shellIsPowerShell)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellPowerShellFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        logical(LK)                     :: shellIsPowerShell
    end function

    impure module function isShellPowerShellFailedMsg(failed, errmsg) result(shellIsPowerShell)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isShellPowerShellFailedMsg
#endif
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        logical(LK)                     :: shellIsPowerShell
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the attempt to run the requested shell command **is successful**.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.
    !>
    !>  \details
    !>  This procedure provides a convenient functional interface to the Fortran intrinsic `execute_command_line()` by
    !>  automatically interpreting the output error codes and reporting the failure as the function output if it occurs.
    !>
    !>  \param[in]      command     :   The input scalar `character` of default kind \SK containing the system shell command to run.
    !>  \param[in]      wait        :   The input scalar `logical` of default kind \LK with the same functionality as the `wait` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  <ul>
    !>                                      <li>    If `.true.`, the procedure will wait for the `command` to terminate before returning the control to the Fortran program.<br>
    !>                                      <li>    Otherwise, the process of running the requested shell command will be done **asynchronously**.<br>
    !>                                  </ul>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[inout]   exitstat    :   The input/output scalar `integer` of default kind \IK with the same functionality as the `exitstat` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  On output, it is assigned the **process exit status** from the command, **only if** `wait = .true.`. If `wait = .false.`, the value of `exitstat` remains untouched.<br>
    !>                                  The value and meaning of this output argument is system dependent, but a non-zero value typically indicates the occurrence of an error within the execution of the specified `command`.<br>
    !>                                  (**optional**. If missing, no exit status will be returned.)
    !>  \param[out]     cmdstat     :   The output scalar `integer` of default kind \IK with the same functionality as the `cmdstat` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  On output, it is assigned<br>
    !>                                  <ul>
    !>                                      <li>     `0` if the process of running the requested shell command terminates without error,
    !>                                      <li>    `-1` if the processor does not support command execution,
    !>                                      <li>    `-2` if `wait = .true.` is specified but the processor does not support asynchronous command execution,
    !>                                      <li>    a positive value if any other error occurs.
    !>                                  </ul>
    !>                                  (**optional**. If missing, no `cmdstat` will be returned.)
    !>  \param[inout]   cmdmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter the same
    !>                                  functionality as the `cmdmsg` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  If present and the output `cmdstat` is assigned a positive value, then `cmdmsg` will be
    !>                                  assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  If the length of `cmdstat` is too short for the output error message, the message tail will be clipped as needed.<br>
    !>                                  A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters for `cmdmsg` is likely enough to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> an error occurs while running the `command`.<br>
    !>                                  The failure is determined by checking for a non-zero value of `cmdstat`.<br>
    !>                                  The specific value of `exitstat` depends on both the processor and the specified value for the input argument `wait`.<br>
    !>                                  Moreover, the processor may execute the specified `command` without failure even though `exitstat` indicates the occurrence of an error within the command.<br>
    !>                                  As such, the output value of `exitstat` has no impact on the value of the output `failed`.<br>
    !>
    !>  \interface{isFailedExec}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedExec
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedExec(command, wait = wait, exitstat = exitstat, cmdstat = cmdstat, cmdmsg = cmdmsg)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure is prone to abuse by requesting to run system commands with devastating consequences.
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>
    !>  \example{isFailedExec}
    !>  \include{lineno} example/pm_sysShell/isFailedExec/main.F90
    !>  \compilef{isFailedExec}
    !>  \output{isFailedExec}
    !>  \include{lineno} example/pm_sysShell/isFailedExec/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \todo
    !>  \plow
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to Fortran intrinsic procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when
    !>  they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedExec}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedExec
    impure elemental module function isFailedExec(command, wait, exitstat, cmdstat, cmdmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedExec
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: command
        logical(LK)     , intent(in)    , optional  :: wait
        integer(IK)     , intent(inout) , optional  :: exitstat
        integer(IK)     , intent(out)   , optional  :: cmdstat
        character(*, SK), intent(inout) , optional  :: cmdmsg
        logical(LK)                                 :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the attempt to fetch
    !>  the shape (height and width) of the current shell terminal is successful.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.<br>
    !>
    !>  \details
    !>  This procedure utilizes shell-specific commands and environment
    !>  variables to get the properties of the current terminal window.
    !>
    !>  \param[in]      width   :   The output scalar `integer` of default kind \IK containing the width (in characters) of the
    !>                              current terminal windows, <b>if and only if</b> the width retrieval is successful.
    !>  \param[in]      height  :   The output scalar `integer` of default kind \IK containing the height (in characters) of the
    !>                              current terminal windows, <b>if and only if</b> the width retrieval is successful.
    !>
    !>  \return
    !>  `failed`                :   The output scalar `logical` of default kind \LK. It is `.true.`
    !>                              <b>if and only if</b> an error occurs while fetching the requested property.<br>
    !>
    !>  \interface{isFailedGetShellShape}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedGetShellShape
    !>      use pm_kind, only: LK, IK
    !>      logical(LK) :: failed
    !>      integer(IK) :: width
    !>
    !>      failed = isFailedGetShellShape(width, height)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [isFailedGetShellHeight](@ref pm_sysShell::isFailedGetShellHeight)<br>
    !>
    !>  \example{isFailedGetShellShape}
    !>  \include{lineno} example/pm_sysShell/isFailedGetShellShape/main.F90
    !>  \compilef{isFailedGetShellShape}
    !>  \output{isFailedGetShellShape}
    !>  \include{lineno} example/pm_sysShell/isFailedGetShellShape/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \todo
    !>  \plow
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to Fortran intrinsic procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when
    !>  they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedGetShellShape}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGetShellShape
    impure module function isFailedGetShellShape(width, height) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetShellShape
#endif
        use pm_kind, only: IKG => IK
        integer(IK)     , intent(out)   :: width
        integer(IK)     , intent(out)   :: height
        logical(LK)                     :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the attempt to fetch the width (in characters) of the current shell terminal is successful.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.
    !>
    !>  \details
    !>  This procedure utilizes shell-specific commands and environment
    !>  variables to get the properties of the current terminal window.
    !>
    !>  \param[in]      width   :   The output scalar `integer` of default kind \IK containing the width (in characters) of the
    !>                              current terminal windows, <b>if and only if</b> the width retrieval is successful.
    !>
    !>  \return
    !>  `failed`                :   The output scalar `logical` of default kind \LK. It is `.true.`
    !>                              <b>if and only if</b> an error occurs while fetching the requested property.<br>
    !>
    !>  \interface{isFailedGetShellWidth}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedGetShellWidth
    !>      use pm_kind, only: LK, IK
    !>      logical(LK) :: failed
    !>      integer(IK) :: width
    !>
    !>      failed = isFailedGetShellWidth(width)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [isFailedGetShellHeight](@ref pm_sysShell::isFailedGetShellHeight)<br>
    !>
    !>  \example{isFailedGetShellWidth}
    !>  \include{lineno} example/pm_sysShell/isFailedGetShellWidth/main.F90
    !>  \compilef{isFailedGetShellWidth}
    !>  \output{isFailedGetShellWidth}
    !>  \include{lineno} example/pm_sysShell/isFailedGetShellWidth/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \todo
    !>  \plow
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to Fortran intrinsic procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when
    !>  they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedGetShellWidth}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGetShellWidth
    impure module function isFailedGetShellWidth(width) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetShellWidth
#endif
        use pm_kind, only: IKG => IK
        integer(IK)     , intent(out)   :: width
        logical(LK)                     :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `failed = .false.` if the attempt to fetch the height (in characters) of the current shell terminal is successful.<br>
    !>  Otherwise, return `failed = .true.` indicating the occurrence of an error.
    !>
    !>  \details
    !>  This procedure utilizes shell-specific commands and environment
    !>  variables to get the properties of the current terminal window.
    !>
    !>  \param[in]      height  :   The output scalar `integer` of default kind \IK containing the height (in characters) of the
    !>                              current terminal windows, <b>if and only if</b> the height retrieval is successful.
    !>
    !>  \return
    !>  `failed`                :   The output scalar `logical` of default kind \LK. It is `.true.`
    !>                              <b>if and only if</b> an error occurs while fetching the requested property.<br>
    !>
    !>  \interface{isFailedGetShellHeight}
    !>  \code{.F90}
    !>
    !>      use pm_sysShell, only: isFailedGetShellHeight
    !>      use pm_kind, only: LK, IK
    !>      logical(LK) :: failed
    !>      integer(IK) :: height
    !>
    !>      failed = isFailedGetShellHeight(height)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [isFailedGetShellWidth](@ref pm_sysShell::isFailedGetShellWidth)<br>
    !>
    !>  \example{isFailedGetShellHeight}
    !>  \include{lineno} example/pm_sysShell/isFailedGetShellHeight/main.F90
    !>  \compilef{isFailedGetShellHeight}
    !>  \output{isFailedGetShellHeight}
    !>  \include{lineno} example/pm_sysShell/isFailedGetShellHeight/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysShell](@ref test_pm_sysShell)
    !>
    !>  \todo
    !>  \plow
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to Fortran intrinsic procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when
    !>  they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedGetShellHeight}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGetShellHeight
    impure module function isFailedGetShellHeight(height) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetShellHeight
#endif
        use pm_kind, only: IKG => IK
        integer(IK)     , intent(out)   :: height
        logical(LK)                     :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sysShell ! LCOV_EXCL_LINE