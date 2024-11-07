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

submodule (pm_sysShell) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_strASCII, only: getStrLower
    use pm_io, only: setContentsFrom
    use pm_io, only: LEN_IOMSG

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `private` module object of type [shell_type](@ref pm_sysShell::shell_type).
    !>
    !>  \details
    !>  This object is set only once throughout the life of the program to avoid costly redundant construction of shell objects.<br>
    !>  It is set by and exclusively used within the routines of this submodule and nowhere else.<br>
    !>  The allocation status of the object is used as an indicator of its initialization.<br>
    !>
    !>  \final{mc_shell}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(shell_type) :: shellinit_type
        logical(LK) :: initialized = .false._LK
    end type
    type(shellinit_type) :: mc_shell

    !>  \brief
    !>  The scalar `private` module object of type `logical` of default kind \LK.
    !>
    !>  \details
    !>  This scalar runtime constant is used to indicate whether the scalar runtime constant [mc_shell](@ref pm_sysShell::mc_shell) is set.<br>
    !>  This indicator could have readily been the allocation status of an `allocatable` mc_shell.<br>
    !>  However, gfortran as of V. 13 appears to have bug not being able to set the components of the object correctly.<br>
    !>
    !>  \final{mc_shellSet}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    !logical(LK) :: mc_shellSet = .false._LK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `allocatable` scalar of type `character` of default kind \SK containing the path to the system temporary directory.<br>
    !>  This runtime module constant is meant to be allocated and set once at runtime and used later throughout the program.<br>
    character(:, SK), save, allocatable :: mc_sysDirTemp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !function isFailedInitShell(errmsg) result(failed)
    !    use pm_kind, only: SKG => SK
    !    character(*,SKG), intent(inout), optional :: errmsg
    !    logical(LK) :: failed
    !    failed = .not. allocated(mc_shell)
    !    if (failed) then
    !        if (present(errmsg)) then
    !            allocate(mc_shell, source = shell_type(failed, errmsg))
    !        else
    !            allocate(mc_shell, source = shell_type(failed))
    !        end if
    !        !mc_shellSet = .not. failed
    !        if (failed .and. allocated(mc_shell)) deallocate(mc_shell)
    !    end if
    !end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedGetDirTemp
        character(LEN_IOMSG, SK) :: errmsg
        errmsg = SK_""
        failed = isFailedGetDirTemp(dirTemp, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedGetDirTempMsg
        use pm_kind, only: SKG => SK
        use pm_sysPath, only: isDir
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedGetDirTemp()"
        integer(IK), parameter :: MAXLEN = max(len(VARENV_DIRTEMP_UNIX), len(VARENV_DIRTEMP_WINDOWS))
        character(MAXLEN, SK), parameter :: VARENV_DIRTEMP(*) = [character(MAXLEN,SKG) :: VARENV_DIRTEMP_UNIX, VARENV_DIRTEMP_WINDOWS]
        integer(IK) :: ienv
#if     OMP_ENABLED
        !$omp critical
#endif
        if (allocated(mc_sysDirTemp)) then
            dirTemp = mc_sysDirTemp
            failed = .false._LK
        else
            ! We cannot infer the runtime shell through `isShellWindows` or similar functions as they create unlimited circular recursive calls.
            ! Instead, we check all environmental variables for all shells, from posix to windows.
            do ienv = 1_IK, size(VARENV_DIRTEMP, 1, IK)
                failed = isFailedGetEnvVar(trim(VARENV_DIRTEMP(ienv)), dirTemp, errmsg) ! returns empty upon error
                failed = failed .or. logical(len_trim(dirTemp, IK) == 0_IK, LK)
                if (.not. failed) exit
            end do
            ! One last try.
            if (failed) dirTemp = SKG_"/tmp"
            failed = .not. isDir(dirTemp, ienv, errmsg)
            failed = failed .or. ienv /= 0_IK
            if (.not. failed) mc_sysDirTemp = dirTemp
        end if
#if     OMP_ENABLED
        !$omp end critical
#endif
    end procedure
    !#if         INTEL_ENABLED
    !#define     FILE_ARG directory = SKG_"/tmp/"
    !#else
    !#define     FILE_ARG file = SKG_"/tmp/"
    !#endif
    !            inquire(FILE_ARG, exist = exists, iostat = i, iomsg = errmsg) ! last forward slash in FILE_ARG is significant.
    !#undef      FILE_ARG
    !            failed = logical(i /= 0_IK, LK) .or. .not. exists
    !            if (failed) then
    !                errmsg = PROCEDURE_NAME//SK_": Failed to fetch shell temporary dir. "//trim(errmsg) ! LCOV_EXCL_LINE
    !                return ! LCOV_EXCL_LINE
    !            end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedGetOutput
        character(LEN_IOMSG, SK) :: errmsg
        errmsg = SK_""
        failed = isFailedGetOutput(command, output, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedGetOutputMsg
        use pm_err, only: getLine
        use pm_kind, only: SKG => SK
        use pm_str, only: isEndedWith
        use pm_arrayStrip, only: getStripped
        use pm_parallelism, only: getImageID
        use pm_distUnif, only: setUnifRand, xoshiro256ssw_type
        type(xoshiro256ssw_type)        :: rng
        character(*, SK), parameter     :: NLC = new_line(SKG_"a"), LF = char(10, SKG), CR = char(13, SKG)
        character(*, SK), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedGetOutputMsg()"
        character(:,SKG), allocatable   :: filetemp, basetemp
        integer(IK)     , parameter     :: lenNLC = len(NLC, IK)
        integer(IK)     , parameter     :: NTRY = 100000_IK
       !integer(IK)                     :: lenBaseTempPlus1
        integer(IK)                     :: lenBaseTemp
        integer(IK)                     :: lenFileTemp
        integer(IK)                     :: counter
       !integer(IK)                     :: length
        integer(IK)                     :: iostat
#define SET_ERRMSG(LINE)\
errmsg = PROCEDURE_NAME//getLine(LINE)//SK_": Failed to fetch command output. "//trim(errmsg)
        ! Determine the shell type.

        !type(shell_type) :: shell
        !shell = shell_type(failed, errmsg)
        !RETURN_IF(failed) ! fpp

        ! Fetch the system temporary dir and create the base of the path.

        failed = isFailedGetDirTemp(basetemp, errmsg)
        if (failed) then
            basetemp = SKG_"isFailedGetOutputMsg.tmp."
        else
        !basetemp = basetemp//shell%dirsep//SKG_"isFailedGetOutputMsg.tmp"
            !   \todo
            !   \warning
            !   The following assumes that all platforms and runtime shells recognize forward slash as a directory separator.
            !   While this is currently the case, it may not be so in the future. A more robust solution may be necessary.
            !   One solution is to remove dependence on the temporary directory and create the file in the current directory.
            basetemp = basetemp//SKG_"/isFailedGetOutputMsg.tmp."
        end if

        ! Generate a unique file to avoid racing conditions in parallel.

        lenBaseTemp = len(basetemp, IK)
        lenFileTemp = lenBaseTemp + 20_IK ! 18 + 12_IK ! date_and_time() is 18 characters.
        if (allocated(filetemp)) deallocate(filetemp) ! \bug gfortran bug in automatic deallocation as of version 11.
        allocate(character(lenFileTemp, SKG) :: filetemp)
        rng = xoshiro256ssw_type(imageID = getImageID())
        filetemp(1 : lenBaseTemp) = basetemp
        do counter = 1_IK, NTRY
            !call setStr(filetemp(lenBaseTempPlus1 : lenFileTemp), length, counter)
            !inquire(file = filetemp(1 : lenBaseTempPlus1 + length - 1), exist = failed, iostat = iostat, iomsg = errmsg)
            !call setUnifRand(filetemp(lenBaseTemp + 19 :), SKG_"a", SKG_"z")
            !call date_and_time(date = filetemp(lenBaseTemp + 1 : lenBaseTemp + 8), time = filetemp(lenBaseTemp + 9 : lenBaseTemp + 18))
            call setUnifRand(rng, filetemp(lenBaseTemp + 1 :), SKG_"a", SKG_"z")
            inquire(file = filetemp, exist = failed, iostat = iostat, iomsg = errmsg)
            failed = failed .or. logical(iostat /= 0_IK, LK)
            if (.not. failed) exit
        end do
        !lenFileTemp = lenBaseTempPlus1 + length - 1_IK

        if (.not. failed) then
            ! Run the command. Beware that the command execution can fail if the command is nonsensical.
            failed = isFailedExec(command//SKG_" 1>"""//filetemp(1 : lenFileTemp)//SKG_""" 2>&1", cmdmsg = errmsg, exitstat = exitstat)
            if (.not. failed) then
                ! Read the command output from file if it exists.
                call setContentsFrom(filetemp(1 : lenFileTemp), contents = output, iostat = iostat, iomsg = errmsg, del = .true._LK)
                failed = logical(iostat /= 0_IK, LK)
                if (failed) then
                    SET_ERRMSG(__LINE__)
                else
                    ! Remove the end-of-file (linefeed) and Carriage Return characters from the captured text,
                    ! because this is artificially added to the command output when written to the file.
                    !if (output(max(1_IK, len(output, IK) - lenNLC + 1_IK) : len(output, IK)) == NLC) output = output(1 : len(output, IK) - lenNLC)
                    !output = output(getSIL(output, new_line(SKG_"a")) : getSIR(output, new_line(SKG_"a")))
                    output = getStripped(output, NLC)
                    loopStrip: do
                        counter = min(1, len(output))
                        if (output(1 : counter) == LF .or. isEndedWith(output, LF) .or. output(1 : counter) == CR .or. isEndedWith(output, CR)) then
                            output = getStripped(getStripped(output, new_line(SKG_"a")), achar(13, SKG)) ! CR must be removed on Windows OS.
                        else
                            exit loopStrip
                        end if
                    end do loopStrip
                end if
            else
                SET_ERRMSG(__LINE__)
            end if
        else
            SET_ERRMSG(__LINE__)
        end if
#undef  SET_ERRMSG
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedExec
        use pm_kind, only: SKG => SK
        use pm_val2str, only: getStr
        character(:, SK), allocatable               :: details
        logical                                     :: wait_def
        integer                                     :: cmdstat_def
        integer                                     :: exitstat_def
        character(*, SK), parameter                 :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedExec()"
        character(*, SK), parameter                 :: LF = new_line(SK_"a")

        if (present(wait)) then
            wait_def = logical(wait, LK)
        else
            wait_def = .true._LK
        end if

        if (present(exitstat)) exitstat_def = int(exitstat)

        if (present(cmdmsg)) then
            call execute_command_line(command, wait_def, exitstat_def, cmdstat_def, cmdmsg)
        else
            call execute_command_line(command, wait_def, exitstat_def, cmdstat_def)
        end if

        if (present(exitstat)) exitstat = int(exitstat_def, IK)
        if (present(cmdstat)) cmdstat = int(cmdstat_def, IK)

        failed = logical(cmdstat_def /= 0, LK)

        if (failed .and. present(cmdmsg)) then
            details = LF//SK_"exitstat = "//getStr(exitstat_def)//LF//SK_"cmdstat = "//getStr(cmdstat_def)//LF//SK_"command = "//getStr(command) ! LCOV_EXCL_LINE
            if (cmdstat_def == -1) then ! LCOV_EXCL_LINE
                cmdmsg = PROCEDURE_NAME//SK_": The processor does not support command execution. command: "//details ! LCOV_EXCL_LINE
            elseif (cmdstat_def == -2 .and. wait_def) then ! LCOV_EXCL_LINE
                cmdmsg = PROCEDURE_NAME//SK_": processor does not support asynchronous command execution. command: "//details ! LCOV_EXCL_LINE
            elseif (cmdstat_def > 0) then ! LCOV_EXCL_LINE
                cmdmsg = PROCEDURE_NAME//SK_": "//trim(adjustl(cmdmsg))//details ! LCOV_EXCL_LINE
            else ! LCOV_EXCL_LINE
                error stop "How on Earth this could happen?! Hey compiler, you are violating the 2018 Fortran standard here." ! LCOV_EXCL_LINE
            end if
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellWindows
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsWindows = isShellWindows(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isShellWindows(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellWindowsFailed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsWindows = isShellWindows(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellWindowsFailedMsg
        type(shell_type) :: shell
        shell = shell_type(failed, errmsg)
        !failed = isFailedInitShell(errmsg)
        if (failed) then
            errmsg = MODULE_NAME//SK_"@isShellWindowsFailedMsg(): Failed to fetch shell type. "//trim(errmsg) ! LCOV_EXCL_LINE
            shellIsWindows = .false._LK
        else
            shellIsWindows = shell%is%windows
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellPosix
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsPosix = isShellPosix(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isShellPosix(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellPosixFailed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsPosix = isShellPosix(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellPosixFailedMsg
        !failed = .not. mc_shellSet
        !if (failed) then
        !    if (allocated(mc_shell)) deallocate(mc_shell); mc_shell = shell_type(failed, errmsg)
        !    mc_shellSet = .not. failed
        !end if
        type(shell_type) :: shell
        shell = shell_type(failed, errmsg)
        !failed = isFailedInitShell(errmsg)
        if (failed) then
            errmsg = MODULE_NAME//SK_"@isShellPosixFailedMsg(): Failed to fetch shell type. "//trim(errmsg) ! LCOV_EXCL_LINE
            shellIsPosix = .false._LK
        else
            shellIsPosix = shell%is%posix
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellCMD
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsCMD = isShellCMD(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isShellCMD(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellCMDFailed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsCMD = isShellCMD(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellCMDFailedMsg
        !failed = .not. mc_shellSet
        !if (failed) then
        !    if (allocated(mc_shell)) deallocate(mc_shell); mc_shell = shell_type(failed, errmsg)
        !    mc_shellSet = .not. failed
        !end if
        !failed = isFailedInitShell(errmsg)
        type(shell_type) :: shell
        shell = shell_type(failed, errmsg)
        if (failed) then
            errmsg = MODULE_NAME//SK_"@isShellCMDFailedMsg(): Failed to fetch shell type. "//trim(errmsg) ! LCOV_EXCL_LINE
            shellIsCMD = .false._LK
        else
            shellIsCMD = shell%is%cmd
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellPowerShell
        logical(LK) :: failed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsPowerShell = isShellPowerShell(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@isShellPowerShell(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellPowerShellFailed
        character(LEN_IOMSG, SK) :: errmsg
        shellIsPowerShell = isShellPowerShell(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isShellPowerShellFailedMsg
        !failed = .not. mc_shellSet
        !if (failed) then
        !    if (allocated(mc_shell)) deallocate(mc_shell); mc_shell = shell_type(failed, errmsg)
        !    mc_shellSet = .not. failed
        !end if
        !failed = isFailedInitShell(errmsg)
        type(shell_type) :: shell
        shell = shell_type(failed, errmsg)
        if (failed) then
            errmsg = MODULE_NAME//SK_"@isShellPowerShellFailedMsg(): Failed to fetch shell type. "//trim(errmsg) ! LCOV_EXCL_LINE
            shellIsPowerShell = .false._LK
        else
            shellIsPowerShell = shell%is%powershell
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure shellis_typer
        character(LEN_IOMSG, SK) :: errmsg
        logical(LK) :: failed
        errmsg = SK_""
        shellis = shellis_type(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@shellis_typer(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure shellis_typerFailed
        character(LEN_IOMSG, SK) :: errmsg
        errmsg = SK_""
        shellis = shellis_type(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure shellis_typerFailedMsg
        !failed = .not. mc_shellSet
        !if (failed) then
        !    if (allocated(mc_shell)) deallocate(mc_shell); mc_shell = shell_type(failed, errmsg)
        !    mc_shellSet = .not. failed
        !end if
        !failed = isFailedInitShell(errmsg)
        type(shell_type) :: shell
        shell = shell_type(failed, errmsg)
        if (failed) then
            errmsg = MODULE_NAME//SK_"@shellis_typerFailedMsg(): "//trim(errmsg) ! LCOV_EXCL_LINE
        else
            shellis = shell%is
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure shell_typer
        character(LEN_IOMSG, SK) :: errmsg
        logical(LK) :: failed
        errmsg = SK_""
        shell = shell_type(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@shell_typer(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure shell_typerFailed
        character(31, SK) :: errmsg
        errmsg = SK_""
        shell = shell_type(failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure shell_typerFailedMsg

        use pm_parallelism, only: getImageID
        use pm_sysPath, only: DIR_SEP_WINDOWS, DIR_SEP_WINDOWS_ALL
        use pm_sysPath, only: PATH_SEP_POSIX, PATH_SEP_WINDOWS
        use pm_sysPath, only: DIR_SEP_POSIX, DIR_SEP_POSIX_ALL
        use pm_distUnif, only: setUnifRand, xoshiro256ssw_type
        use pm_arrayStrip, only: getSIL, getSIR
        use pm_kind, only: SKG => SK
        use pm_err, only: getLine

        type(xoshiro256ssw_type)        :: rng
        character(*, SK), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@shell_typer()"
        character(*, SK), parameter     :: LF = new_line(SK_"a")
        character(:, SK), allocatable   :: filetemp, basetemp
        integer(IK)     , parameter     :: NTRY = 100000_IK
        integer(IK)                     :: iostat, counter
        integer(IK)                     :: lenBaseTemp
        integer(IK)                     :: lenFileTemp
#define SET_ERRMSG(LINE)\
errmsg = PROCEDURE_NAME//getLine(LINE)//SK_": Failed to fetch shell name. "//trim(errmsg)

        if (mc_shell%initialized) then

            shell = mc_shell%shell_type
            failed = .false._LK

        else

            ! Initialize the allocatable components to the default.

            shell%name = SKG_""
            shell%dirseps = DIR_SEP_POSIX_ALL

            ! Fetch the system temporary dir and create the base of the path.

            failed = isFailedGetDirTemp(basetemp, errmsg)
            if (failed) then
                basetemp = SKG_"shell_typer.tmp."
            else
                !   \todo
                !   \warning
                !   The following assumes that all platforms and runtime shells recognize forward slash as a directory separator.
                !   While this is currently the case, it may not be so in the future. A more robust solution may be necessary.
                !   One solution is to remove dependence on the temporary directory and create the file in the current directory.
                basetemp = basetemp//SKG_"/shell_typer.tmp."
            end if

            ! Generate a unique file name to avoid racing conditions in parallel.

            lenBaseTemp = len(basetemp, IK)
            lenFileTemp = lenBaseTemp + 21 ! + 12 ! + 1 dot + date and time length is 18.
            if (allocated(filetemp)) deallocate(filetemp) ! \bug gfortran bug in automatic deallocation as of version 11.
            allocate(character(lenFileTemp, SKG) :: filetemp)
            rng = xoshiro256ssw_type(imageID = getImageID())
            filetemp(1 : lenBaseTemp) = basetemp
            do counter = 1_IK, NTRY
                call setUnifRand(rng, filetemp(lenBaseTemp + 1 :), SKG_"a", SKG_"z")
                inquire(file = filetemp, exist = failed, iostat = iostat, iomsg = errmsg)
                failed = failed .or. logical(iostat /= 0_IK, LK)
                if (.not. failed) exit
            end do

            if (failed) then
                SET_ERRMSG(__LINE__)
            else
                ! First check for Unix shells, as the command does not lead to odd syntax errors on a Windows shell.
                ! $SHELL gives the full path to the default shell (not the runtime shell).
                ! $0 gives the name (OR the path, e.g., in Git Bash) of the current shell.
                ! $0 works in the following shells as of 2022: ash, bash, csh, dash, sh, tcsh, zsh, yash
                ! $0 does not work in the following shells as of 2022: cmd, fish, powershell, pwsh,
                ! WARNING
                ! You may think `isFailedExec()` can be replaced with `isFailedGetOutput()` in the following.
                ! That is a false premise.
                ! Such a replacement creates cyclic dependence on shell_type constructor which then fails. Do not replace it.
                failed = isFailedExec(SK_"echo $0 1>"""//filetemp//""" 2>&1")
                !write(*,*) filetemp
                call setContentsFrom(file = filetemp, contents = shell%name, iostat = iostat, iomsg = errmsg)!, del = .true._LK)
                failed = logical(iostat /= 0_IK, LK)
                if (failed) then
                    SET_ERRMSG(__LINE__) ! LCOV_EXCL_LINE
                else
                    if (0 < len(shell%name)) shell%name = shell%name(getSIL(shell%name, LF) : getSIR(shell%name, LF))
                    ! Assign the component values based on the contents of the command output.
                    call setShell()
                    if (shell%is%posix .and. shell%is%sh) then
                        ! If shell is sh, it may be a symlink, or some other shell renamed. If so, get the target shell type.
                        failed = isFailedExec(SK_"command -v sh 1>"""//filetemp//""" 2>&1")
                        if (.not. failed) then
                            call setContentsFrom(file = filetemp, contents = shell%name, iostat = iostat, iomsg = errmsg, del = .true._LK)
                            shell%name = shell%name(getSIL(shell%name, LF) : getSIR(shell%name, LF))
                            failed = logical(iostat /= 0_IK, LK)
                        end if
                        ! Fetch the symlink target.
                        if (.not. failed) then
                            failed = isFailedExec(SKG_"ls -l """//shell%name//SKG_""" 1>"""//filetemp//""" 2>&1")
                            if (.not. failed) then
                                call setContentsFrom(file = filetemp, contents = shell%name, iostat = iostat, iomsg = errmsg, del = .true._LK)
                                failed = logical(iostat /= 0_IK, LK)
                                if (.not. failed) call setShell()
                            end if
                        end if
                        if (failed) then
                            errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg) ! LCOV_EXCL_LINE
                            shell%name = SKG_"sh" ! LCOV_EXCL_LINE
                            failed = .false._LK ! LCOV_EXCL_LINE
                        end if
                    elseif (.not. shell%is%posix) then
                        ! First, check for fish shell.
                        failed = isFailedGetEnvVar("FISH_VERSION", shell%name, errmsg, length = 63_IK)
                        if (failed) then
                            errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg) ! LCOV_EXCL_LINE
                            shell%name = SKG_"" ! LCOV_EXCL_LINE
                        elseif (len_trim(shell%name, IK) > 0_IK) then
                            !shell%is%posix = .true._LK ! fish is not posix-compliant. Dont you dare setting this to .true.!
                            shell%is%fish = .true._LK
                            shell%name = SKG_"fish"
                        end if
                    end if
                    !  Test for Windows-based terminals.
                    if (.not. failed) then
                        if (.not. (shell%is%posix .or. shell%is%fish)) then
                            ! \note
                            ! Note PowerShell Core pwsh.exe is also available on Unix systems and recognizes the same syntax as in Windows PowerShell BUT NOT ALL.
                            ! It also recognizes POXIS commands but not all!
                            ! Even worse, not all Windows PowerShell commands are defined in Unix PowerShell Core.
                            ! This is super-confusing. What the hell Microsoft are you doing?!
                            ! Fortunately the command for CMD shell recognition is platform-independent.
                            ! LCOV_EXCL_START
                            failed = isFailedExec(SK_"(dir 2>&1 *`|echo CMD >"""//filetemp//""");&<# rem #>echo PowerShell >"""//filetemp//""" 2>&1")
                            call setContentsFrom(file = filetemp, contents = shell%name, iostat = iostat, iomsg = errmsg, del = .true._LK)
                            failed = logical(iostat /= 0_IK, LK)
                            if (failed) SET_ERRMSG(__LINE__) ! LCOV_EXCL_LINE
                            if (.not. failed) then
                                if (index(shell%name, SKG_"CMD", kind = IK) > 0_IK) then
                                    shell%is%windows = .true._LK
                                    shell%is%cmd = .true._LK
                                elseif (index(shell%name, SKG_"PowerShell", kind = IK) > 0_IK) then
                                    shell%is%powershell = .true._LK
#if                                 WINDOWS_ENABLED
                                    shell%is%windows = .true._LK
#elif                               DARWIN_ENABLED || LINUX_ENABLED
                                    shell%is%posix = .true._LK
#else
                                    failed = isFailedGetEnvVar(SK_"OS", basetemp, errmsg, length = 10_IK)
                                    if (failed) then
                                        SET_ERRMSG(__LINE__)
                                    else
                                        shell%is%windows = index(getStrLower(basetemp), SK_"windows") /= 0)
                                        shell%is%posix = .not. shell%is%windows
                                    end if
#endif
                                end if
                            end if
                        end if
                        if (.not. failed) then
                            if (shell%is%windows) then
                                shell%dirseps = DIR_SEP_WINDOWS_ALL
                                shell%pathsep = PATH_SEP_WINDOWS
                                shell%dirsep = DIR_SEP_WINDOWS
                            elseif (shell%is%posix .or. shell%is%fish) then
                                shell%dirseps = DIR_SEP_POSIX_ALL
                                shell%pathsep = PATH_SEP_POSIX
                                shell%dirsep = DIR_SEP_POSIX
                            else
                                SET_ERRMSG(__LINE__) ! LCOV_EXCL_LINE
                                failed = .true._LK ! LCOV_EXCL_LINE
                            end if
                            ! Cache the results.
                            mc_shell%shell_type = shell
                        end if
                    end if
                end if
            end if
        end if

    contains

        subroutine setShell()
            ! Do **not** change the order of the following name checks.
            if (index(shell%name, SKG_"bash", kind = IK) > 0_IK) then
                shell%name = SKG_"bash"
                shell%is%sh = .false._LK
                shell%is%bash = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"dash", kind = IK) > 0_IK) then
                shell%name = SKG_"dash"
                shell%is%sh = .false._LK
                shell%is%dash = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"yash", kind = IK) > 0_IK) then
                shell%name = SKG_"yash"
                shell%is%sh = .false._LK
                shell%is%yash = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"ash", kind = IK) > 0_IK) then
                shell%name = SKG_"ash"
                shell%is%sh = .false._LK
                shell%is%ash = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"tcsh", kind = IK) > 0_IK) then
                shell%name = SKG_"tcsh"
                shell%is%sh = .false._LK
                shell%is%tcsh = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"csh", kind = IK) > 0_IK) then
                shell%name = SKG_"csh"
                shell%is%sh = .false._LK
                shell%is%csh = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"ksh", kind = IK) > 0_IK) then
                shell%name = SKG_"ksh"
                shell%is%sh = .false._LK
                shell%is%ksh = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"zsh", kind = IK) > 0_IK) then
                shell%name = SKG_"zsh"
                shell%is%sh = .false._LK
                shell%is%zsh = .true._LK
                shell%is%posix = .true._LK
            elseif (index(shell%name, SKG_"sh", kind = IK) > 0_IK) then
                shell%name = SKG_"sh"
                shell%is%sh = .true._LK
                shell%is%posix = .true._LK
            end if
        end subroutine

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedGetEnvVar
        character(LEN_IOMSG, SK) :: errmsg
        errmsg = SK_""
        failed = isFailedGetEnvVarMsg(name, value, errmsg, length)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedGetEnvVarMsg

        use pm_kind, only: SKG => SK
        use pm_val2str, only: getStr
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedGetEnvVarMsg()"
        integer :: status, lengthReal
        integer(IK) :: lenVal

        if (present(length)) then
            CHECK_ASSERTION(__LINE__, 0_IK < length, PROCEDURE_NAME//SK_": The input `length` value must be a positive integer. length = "//getStr(length)) ! fpp
            lenVal = length
        else
            lenVal = 8191_IK ! 2**13 - 1
        end if

        allocate(character(lenVal, SKG) :: value, stat = status)
        failed = logical(status /= 0, LK)
        if (failed) then
            errmsg = PROCEDURE_NAME//SK_": Fortran runtime error: Allocation of value failed with size = "//getStr(lenVal)//SK_". stat = "//getStr(status)//SK_"." ! LCOV_EXCL_LINE
            !value = SKG_"" ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        failed = logical(name == SKG_"", LK)
        if (failed) then
            errmsg = PROCEDURE_NAME//SK_": Fortran runtime error: Zero-length string passed as name to get_environment_variable()."
            !value = SK_""
            return
        end if

        loopAdjustLength: do
            call get_environment_variable(name, value, lengthReal, status)
            if (status == 0) then
                if (lengthReal < len(value)) value = value(1 : lengthReal)
                return
            elseif (status == +1) then ! the environment variable does not exist.
                failed = .false._LK
                !value = SK_""
                return
            elseif (status == -1) then ! the value argument is present but too short.
                deallocate(value)
                lenVal = lenVal * 2_IK
                allocate(character(lenVal, SK) :: value)
                cycle loopAdjustLength
            elseif (status == +2) then
                failed = .true._LK ! LCOV_EXCL_LINE
                errmsg = PROCEDURE_NAME//SK_": Failed to fetch the value of environment variable """//name//""". The processor does not support environment variables. status = "//getStr(status) ! LCOV_EXCL_LINE
                !value = SK_"" ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            elseif (status > +2) then
                failed = .true._LK ! LCOV_EXCL_LINE
                errmsg = PROCEDURE_NAME//SK_": Unknown error occurred while fetching the value of the environment variable """//name//""". status = "//getStr(status) ! LCOV_EXCL_LINE
                !value = SK_"" ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            else
                error stop "How on Earth such an error value could happen?! Hey compiler, you have violated the 2018 Fortran Standard. status = "//getStr(status) ! LCOV_EXCL_LINE
            end if
        end do loopAdjustLength

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedPutEnvVar
        ! \warning : The input `dir` must be of default character kind.
        use iso_c_binding, only: c_null_char, c_int
#if     WINDOWS_ENABLED
        interface ! See WIN32 API for details of `_putenv_s`, `_wputenv_s`.
            function setEnvVar(name, value) result(stat) bind(C, name = "_putenv_s")
                use iso_c_binding, only: c_char, c_int
                character(1,c_char), intent(in) :: name(*), value(*)
                integer(c_int) :: stat
            end function
        end interface
        failed = logical(setEnvVar(name//c_null_char, value//c_null_char) /= int(0, c_int), LK)
#else
        integer(c_int), parameter :: OVERWRITE_ENABLED = int(1, c_int)
        interface
            function setEnvVar(name, value, overwrite) result(stat) bind(C, name = "setenv")
                use iso_c_binding, only: c_char, c_int
                character(1,c_char), intent(in) :: name(*), value(*)
                integer(c_int), intent(in) :: overwrite
                integer(c_int) :: stat
            end function
        end interface
        failed = logical(setEnvVar(name//c_null_char, value//c_null_char, OVERWRITE_ENABLED) /= int(0, c_int), LK)
#endif
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedGetShellShape_ENABLED 1
    module procedure isFailedGetShellShape
#include "pm_sysShell@routines.inc.F90"
    end procedure
#undef isFailedGetShellShape_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedGetShellWidth_ENABLED 1
    module procedure isFailedGetShellWidth
#include "pm_sysShell@routines.inc.F90"
    end procedure
#undef isFailedGetShellWidth_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedGetShellHeight_ENABLED 1
    module procedure isFailedGetShellHeight
#include "pm_sysShell@routines.inc.F90"
    end procedure
#undef isFailedGetShellHeight_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
