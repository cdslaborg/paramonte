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
!>  This file contains procedure implementations of [pm_sysPath](@ref pm_sysPath).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sysPath) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define ERROR_STOP_IF(FAILED,MSG) if (FAILED) error stop MSG
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define ERROR_STOP_IF(FAILED,MSG)
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use, intrinsic :: iso_c_binding, only: c_int32_t, c_char, c_null_char
    use pm_os, only: isWindows, isLinux, isDarwin
    use pm_sysShell, only: shellis_type, shell_type
    use pm_sysShell, only: isFailedGetOutput
    use pm_sysShell, only: isFailedGetEnvVar
    use pm_arrayReplace, only: getReplaced
    use pm_arrayReplace, only: setReplaced
    use pm_arrayInsert, only: setInserted
    use pm_arrayRemove, only: getRemoved
    use pm_arrayResize, only: setResized
    use pm_sysShell, only: isFailedExec
    use pm_strASCII, only: isCharAlpha
    use pm_strASCII, only: setStrLower
    use pm_strASCII, only: getStrLower
    use pm_arraySplit, only: setSplit
    use pm_container, only: css_type
    use pm_str, only: isEndedWith
    use pm_str, only: getCharVec
    use pm_val2str, only: getStr

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module procedure constructPathList
!        use pm_kind, only: SKC => SK
!        use pm_val2str, only: getStr
!        use pm_sysShell, only: shellis_type
!        character(*, SK), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@constructPathList()"
!        character(:, SK), allocatable   :: command
!        !logical(LK) :: namesorted, datesorted, sizesorted
!
!        type(shellis_type) :: shellis
!        if (present(failed) .and. present(errmsg)) then
!            shellis = shellis_type(failed, errmsg)
!            errmsg = MODULE_NAME//SK_"@constructPathList(): "//trim(errmsg)
!            if (failed) return ! LCOV_EXCL_LINE
!        elseif (present(failed)) then
!            shellis = shellis_type(failed)
!            if (failed) return ! LCOV_EXCL_LINE
!        else
!            shellis = shellis_type()
!        end if
!
!        ! Determine sorting method.
!
!        if (present(sort)) then
!            if (sort == SK_"name") then
!                namesorted = .true._LK
!                datesorted = .false._LK
!                sizesorted = .false._LK
!            elseif (sort == SK_"date") then
!                namesorted = .false._LK
!                datesorted = .true._LK
!                sizesorted = .false._LK
!            elseif (sort == SK_"date") then
!                namesorted = .false._LK
!                datesorted = .true._LK
!                sizesorted = .false._LK
!            else
!                if (present(failed)) then
!                    failed = .true._LK
!                    if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Unrecognized value for the `sort` input argument: "//getStr(sort)
!                else
!                    error stop PROCEDURE_NAME//SK_": Unrecognized value for the `sort` input argument: "//getStr(sort)
!                end if
!            end if
!        end if
!
!        ! Determine item types to list
!
!        ! Generate the list of items
!
!        !   `ls` flags:
!        !   -a, --all                   :   do not ignore entries starting with `.`.
!        !   -Q, --quote-name            :   enclose entry names in double quotes.
!        !   -A, --almost-all            :   do not list implied `.` and `..`.
!        !   -b, --escape                :   print C-style escapes for nongraphic characters.
!        !   -c with -lt                 :   sort by, and show, ctime (time of last modification of file status information); with -l: show ctime and sort by name; otherwise: sort by ctime, newest first.
!        !   -F, --classify              :   append indicator (one of */=>@|) to entries. Display
!        !                                       -#  a slash (/) immediately after each pathname that is a directory,
!        !                                       -#  an asterisk (`*') after each that is executable,
!        !                                       -#  an at sign (`@') after each symbolic link,
!        !                                       -#  a percent sign (`%') after each whiteout,
!        !                                       -#  an equal sign (`=') after each socket,
!        !                                       -#  avertical bar (`|') after each that is a FIFO.
!        !   --group-directories-first   :   group directories before files; can be augmented with a `--sort` option, but any use of `--sort=none` (`-U`) disables grouping.
!        !   -h, --human-readable        :   with -l and -s, print sizes like 1K 234M 2G etc.
!        !   -k, --kibibytes             :   default to 1024-byte blocks for disk usage; used only with `-s` and per directory totals.
!        !   -m                          :   fill width with a comma separated list of entries.
!        !   -Q, --quote-name            :   enclose entry names in double quotes.
!        !   --quoting-style=WORD        :   use quoting style WORD for entry names: `literal`, `locale`, `shell`, `shell-always`, `shell-escape`, `shell-escape-always`, `c`, `escape` (overrides QUOTING_STYLE environment variable)
!        !   -r, --reverse               :   reverse order while sorting
!
!
!        if (shellis%posix .or. shellis%fish) then
!            if (namesorted)
!            command = "ls -b"
!        elseif (shellis%windows) then
!        else
!            if (present(failed)) then
!                failed = .true._LK
!                if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Unrecognized runtime shell."
!            else
!                error stop PROCEDURE_NAME//SK_": Unrecognized runtime shell."
!            end if
!        end if
!
!    end procedure
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathMatch_ENABLED 1

    module procedure getPathMatch
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef getPathMatch_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathMatch_ENABLED 1

    module procedure setPathMatch
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef setPathMatch_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isDir_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    module procedure isDirDD
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    module procedure isDirII
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isDir_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFile_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    module procedure isFileDD
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    module procedure isFileII
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isFile_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isExtant_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    module procedure isExtantDD
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define II_ENABLED 1

    module procedure isExtantII
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef II_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isExtant_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathVerbatim

        use pm_kind, only: SKC => SK
        use pm_sysShell, only: shellis_type

        type(shellis_type) :: shellis
        if (present(failed) .and. present(errmsg)) then
            shellis = shellis_type(failed, errmsg)
            errmsg = MODULE_NAME//SK_"@getPathVerbatim(): "//trim(errmsg)
            if (failed) return ! LCOV_EXCL_LINE
        elseif (present(failed)) then
            shellis = shellis_type(failed)
            if (failed) return ! LCOV_EXCL_LINE
        else
            shellis = shellis_type()
        end if

        ! The order of the following conditions must be preserved.
        !print *, shellis
        if (shellis%powershell) then
            pathVerbatim = getPathVerbatimPowerShell(path)
        elseif (shellis%posix) then
            pathVerbatim = getPathVerbatimPosix(path)
        elseif (shellis%fish) then
            pathVerbatim = getPathVerbatimFish(path) ! LCOV_EXCL_LINE
        elseif (shellis%cmd) then ! .and. shellis%windows) then
            pathVerbatim = getPathVerbatimCMD(path) ! LCOV_EXCL_LINE
        elseif (present(failed)) then
            failed = .true._LK ! LCOV_EXCL_LINE
            if (present(errmsg)) errmsg = MODULE_NAME//SK_"@getPathVerbatim(): Failed to infer the runtime shell type." ! LCOV_EXCL_LINE
        else
            error stop MODULE_NAME//SK_"@getPathVerbatim(): Failed to infer the runtime shell type." ! LCOV_EXCL_LINE
        end if

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getIndexDirName_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getIndexDirNameDef
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getIndexDirNamePM
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getIndexDirName_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDirName_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirNameDef
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirNamePM
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDirName_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getIndexBaseName_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getIndexBaseNameDef
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getIndexBaseNamePM
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getIndexBaseName_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBaseName_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getBaseNameDef
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getBaseNamePM
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBaseName_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module procedure getIndexBaseName
!        integer(IK) :: offset
!        GET_OFFSET ! fpp
!        indexBaseName = scan(path(1 : offset-1), dirsep, back = .true., kind = IK) + 1_IK
!    end procedure
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module procedure getIndexBaseNamePM
!        indexBaseName = scan(path(1:len(path,IK)), dirsep, back = .true., kind = IK) + 1_IK
!    end procedure
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module procedure getBaseName
!        use pm_kind, only: SKC => SK
!        integer(IK) :: offset
!        integer(IK) :: lenPath
!        lenPath = len(path, IK)
!        if (lenPath > 1_IK) then
!            GET_OFFSET ! fpp
!            basename = path(getIndexBaseName(path(1 : offset), dirsep): offset)
!        else
!            basename = path
!        end if
!    end procedure
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module procedure getBaseNamePM
!        if (asis) then
!            basename = path(getIndexBaseName(path, dirsep, asis):)
!        else
!            basename = getBaseName(path, dirsep)
!        end if
!    end procedure
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getExtName
        extname = path(getIndexExtName(path, dirsep):)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getFileName
        filename = path(getIndexBaseName(path, dirsep) : getIndexExtName(path, dirsep) - 1_IK)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getIndexExtName
        use pm_kind, only: SKC => SK
        indexExtName = len(path, IK)
        if (indexExtName > 0_IK) then
            loopSearch: do
                if (path(indexExtName:indexExtName) == SKC_".") then
                    return
                elseif (path(indexExtName:indexExtName) == dirsep .or. indexExtName == 1_IK) then
                    indexExtName = len(path, IK) + 1_IK
                    return
                end if
                indexExtName = indexExtName - 1_IK
            end do loopSearch
        else
            indexExtName = indexExtName + 1_IK
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
    module procedure getPathDensePosix
        use pm_kind, only: SKC => SK
        character(*,SKC), parameter     :: PSP = DIR_SEP_POSIX
        character(*,SKC), parameter     :: CDIR = PSP//"."//PSP ! current dir pattern.
        character(*,SKC), parameter     :: PDIR = PSP//".."//PSP ! parent dir pattern.
        integer(IK)                     :: istart, iend, j, lenPath
        lenPath = len(path, IK)
        allocate(character(lenPath,SKC) :: pathDense)
        CHECK_ASSERTION(__LINE__, len(PSP, IK) == 1_IK, SK_"@getPathDensePosix(): Posix directory separator must be a single character.")
        istart = 1_IK
        iend = istart
        do
            if (iend < lenPath - 2_IK) then
                if (path(i:i+1) == CDIR) then
                    i = i + 2_IK
                elseif (path(i:i+1_IK) == "/../"//PSP) then
                    i = i + 3_IK
                end if
            else
                exit
            end if
        end do
        if (lenPath > 3_IK) then
        elseif (lenPath == 3_IK) then
            if (
        elseif (lenPath == 2_IK) then
        elseif (lenPath == 1_IK) then
        end if
    end procedure
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathAbs
        logical(LK) :: failed
        character(2047, SK) :: errmsg
        errmsg = SK_""
        pathAbs = getPathAbs(path, failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@getPathAbs(): "//trim(errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathAbsFailed
        character(2047, SK) :: errmsg
        errmsg = SK_""
        pathAbs = getPathAbs(path, failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathAbsFailedMsg

        use pm_sysShell, only: isFailedGetOutput
        use pm_sysShell, only: shellis_type
        use pm_kind, only: SKC => SK
        type(shellis_type) :: shellis

        failed = .false._LK
        blockNonEmptyPath: if (len(path,IK) > 0_IK) then
!#if         INTEL_ENABLED
!            !   Try the intel compiler approach.
!            !   \warning
!            !   The extension `FullPathQQ()` is apparently sensitive to last dirsep,
!            !   causing a doubling of the basename of the input `path` if it ends with a dirsep.
!            !   For now, the bypass is to avoid the use of this Intel compiler extension as its behavior seems ambiguous.
!            block
!                use ifport, only: FullPathQQ
!                integer(IK) :: lenPath
!                allocate(character(MAX_LEN_FILE_PATH,SKC) :: pathAbs)
!                lenPath = FullPathQQ(path, pathAbs)
!                if (0_IK < lenPath) then
!                    pathAbs = pathAbs(1 : lenPath)
!                    return
!                end if
!                !failed = logical(lenPath > 0_IK, LK)
!                !if (failed) errmsg = MODULE_NAME//SK_"@getPathAbsFailedMsg(): Failed to resolve the input path." ! LCOV_EXCL_LINE
!                ! try other methods
!            end block
!#endif
            shellis = shellis_type(failed, errmsg)

            ! Try the Microsoft PowerShell approach.

            if (shellis%powershell .or. shellis%cmd) then
                ! LCOV_EXCL_START
                block

                    !logical(LK) :: endsWithNLC
                    !character(1,SKC) :: APPENDIX = SKC_" "
                    character(:,SKC), allocatable :: command
                    character(*,SKC), parameter :: NLC = new_line(SKC_"a")

                    !   \warning
                    !   The presence of trailing new line character in the input path can be problematic when the absolute path is inferred via Microsoft PowerShell.<br>
                    !   By default, the procedures of this generic interface strip any leading or trailing new line characters via the call to [isFailedGetOutput](@ref pm_sysShell::isFailedGetOutput).<br>
                    !   This issue becomes relevant only if Microsoft PowerShell Core is the default runtime shell on a POSIX-compliant shell, which is as of 2022 almost never the case.
                    !   Update: The above behavior of `isFailedGetOutput()` was a bug that is now resolved. As such, this discussion is irrelevant.

                    !endsWithNLC = logical(path(len(path, IK) - len(NLC, IK) + 1_IK :) == NLC, LK)
                    !if (endsWithNLC) APPENDIX = SKC_"X" ! ensure the significant new line character is not stripped incorrectly by `isFailedGetOutput()`.

                    !   The order of the conditionals in the following should not change (powershell must be first).

                    if (shellis%powershell) then
                        !command = SKC_"$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath('"//path//trim(APPENDIX)//SKC_"')" ! LCOV_EXCL_LINE
                        command = SKC_"$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath('"//path//SKC_"')" ! LCOV_EXCL_LINE
                    elseif (shellis%cmd) then
                        !command = SKC_"powershell -command ""$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath('"//path//trim(APPENDIX)//SKC_"')""" ! LCOV_EXCL_LINE
                        command = SKC_"powershell -command ""$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath('"//path//SKC_"')""" ! LCOV_EXCL_LINE
                    end if

                    failed = isFailedGetOutput(command, pathAbs, errmsg)
                    !if (failed) errmsg = MODULE_NAME//SK_"@getPathAbsFailedMsg(): Failed to fetch powershell command output. "//trim(errmsg) ! LCOV_EXCL_LINE
                    if (.not. failed) return
                    !if (.not. failed) then
                    !    if (endsWithNLC) pathAbs = pathAbs(1 : len(path, IK) - 1)
                    !    return
                    !end if

                end block
                ! LCOV_EXCL_STOP
            end if

            ! Try the POSIX realpath() approach.
            ! The realpath() approach is too complex and involved and its result is not compatible with the Windows approach (that does not resolve symlinks).
            !#if     !WINDOWS_ENABLED
            !! The Unix fencing is required because `realpath` does not exist on Windows OS.
            !elseif (shellis%posix) then
            !    block
            !        use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_null_char, c_associated
            !        character(*,SKC), parameter     :: PSP = DIR_SEP_POSIX
            !        character(:,SKC), allocatable   :: relPathStr, absPathStr
            !        integer(IK)     , parameter     :: NTRY = 8_IK ! max number of times extending the buffer.
            !        integer(IK)                     :: absPathLen, maxPathLen, itry
            !        type(c_ptr)                     :: ptr
            !        interface
            !            function realpath(path,resolved_path) result(absPathStrPtr) bind(c, name = "realpath")
            !                use, intrinsic :: iso_c_binding
            !                character(1, c_char), intent(in)    :: path(*)
            !                character(1, c_char), intent(out)   :: resolved_path(*)
            !                type(c_ptr)                         :: absPathStrPtr
            !            end function realpath
            !        end interface
            !        ERROR _STOP _IF(DIR_SEP_POSIX /= DIR_SEP_POSIX_ALL, SK_"@isPathAbsPosix(): Internal error occurred. `DIR_SEP_POSIX == DIR_SEP_POSIX_ALL` must hold."
            !        if (isPathAbsPosix(path)) then; pathAbs = path; return; end if
            !        relPathStr = getPathDense(path, PSP) ! no compacting should be done as it can destroy symlinks.
            !        loopFindExistingDirName: do
            !            maxPathLen = 1024_IK
            !            loopSufficientPathLen: do itry = 1_IK, NTRY
            !                absPathStr = repeat(SKC_" ", maxPathLen)
            !                absPathStrPtr = realpath(relPathStr//c_null_char, absPathStr)
            !                ! Determine the first null char
            !                do absPathLen = 1_IK, maxPathLen
            !                   if(absPathStr(absPathLen:absPathLen) /= c_null_char) cycle
            !                   absPathLen = absPathLen - 1_IK
            !                   exit
            !                end do
            !                if (absPathLen < maxPathLen) exit loopSufficientPathLen
            !                maxPathLen = maxPathLen * 2_IK
            !            end do loopSufficientPathLen
            !            failed = itry > NTRY
            !            if (failed) then
            !                errmsg = MODULE_NAME//SK_"@getPathAbsFailedMsg(): Failed to detect the null character in the output of realpath(). This is highly unusual." ! LCOV_EXCL_LINE
            !                ptr = c_null_ptr ! LCOV_EXCL_LINE
            !                return ! LCOV_EXCL_LINE
            !            end if
            !            if (c_associated(absPathStrPtr)) exit loopFindExistingDirName
            !        end do loopFindExistingDirName
            !        ptr = c_null_ptr
            !    end block
            !#endif

            ! Ensure the input non-empty path is indeed relative.

            if (isPathAbsPosix(path) .or. (shellis%windows .and. isPathAbsWindows(path))) then
                pathAbs = path
                return
            end if

        end if blockNonEmptyPath

        ! We get here only if the input path is empty or all other approaches have failed.

        block
            character(:,SKC), allocatable   :: dirc
            character(1,SKC), parameter     :: PSP = DIR_SEP_POSIX
            character(1,SKC), parameter     :: PSW = DIR_SEP_WINDOWS
            dirc = getDirCurrent(failed, errmsg)
            if (.not. failed) then
                if (len(path, IK) == 0_IK) then
                    pathAbs = dirc
                elseif (shellis%windows) then
                    pathAbs = dirc//PSW//path
                else
                    pathAbs = dirc//PSP//path
                end if
                return
            end if
        end block

        ERROR_STOP_IF(.not. failed, MODULE_NAME//SK_"@getPathAbsFailedMsg(): Internal library error: The program should have failed by this point.")
        errmsg = MODULE_NAME//SK_"@getPathAbsFailedMsg(): "//trim(errmsg) ! LCOV_EXCL_LINE
        pathAbs = path

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirCurrent
        character(2047, SK) :: errmsg
        logical(LK) :: failed
        dirCurrent = getDirCurrent(failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@getDirCurrent(): "//trim(errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirCurrentFailed
        use pm_kind, only: SKC => SK
#if     INTEL_ENABLED || GNU_ENABLED
#if     INTEL_ENABLED
        use ifport, only: getcwd
#endif
        use pm_arrayResize, only: setResized
        integer(IK) :: dirCurrentLen
        dirCurrentLen = 127_IK
        do
            call setResized(dirCurrent, dirCurrentLen)
            ! This yields an ICE with gfortran on macOS arm64 M2.
            ! It seems like external routines `getcwd()` cannot be passed directly to intrinsic functions.
            ! failed = logical(getcwd(dirCurrent) /= 0, LK)
            failed = getcwd(dirCurrent) /= 0
            if (.not. failed) exit
            if (MAX_LEN_FILE_PATH < dirCurrentLen) return
            dirCurrentLen = 2_IK * (dirCurrentLen + 1_IK)
            cycle
        end do
        if (failed) then
            dirCurrent = SKC_"."
        else
            dirCurrent = trim(dirCurrent)
        end if
#else
        character(2047, SK) :: errmsg
        dirCurrent = getDirCurrent(failed, errmsg)
#endif
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirCurrentFailedMsg
        ! the environment variables CD and perhaps PWD are unreliable and may not be defined. Use the command form instead.
        use pm_sysShell, only: isFailedGetOutput
        use pm_sysShell, only: shellis_type
        use pm_kind, only: SKC => SK
        type(shellis_type) :: shellis
        shellis = shellis_type(failed, errmsg)
        if (shellis%powershell) then
            failed = isFailedGetOutput(SKC_"echo (pwd).path", dirCurrent, errmsg)
        elseif (shellis%posix .or. shellis%fish) then
            failed = isFailedGetOutput(SKC_"pwd", dirCurrent, errmsg)
        elseif (shellis%cmd) then
            failed = isFailedGetOutput(SKC_"cd", dirCurrent, errmsg)
        end if
        if (failed) then
            errmsg = MODULE_NAME//SK_"@getDirCurrent(): Failed to fetch `pwd` command output. "//trim(errmsg) ! LCOV_EXCL_LINE
            dirCurrent = SKC_"." ! LCOV_EXCL_LINE
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirHome
        use pm_kind, only: SKC => SK
        use pm_val2str, only: getStr
        use pm_sysShell, only: shellis_type
        use pm_sysShell, only: isFailedGetEnvVar, isFailedGetOutput
        character(:,SKC), allocatable   :: homeDrive, homePath, user_def
        character(*,SKC), parameter     :: PSWA = DIR_SEP_WINDOWS_ALL
        character(*,SKC), parameter     :: PSPA = DIR_SEP_POSIX_ALL
        character(*,SKC), parameter     :: PSW = DIR_SEP_WINDOWS
        character(*,SKC), parameter     :: PSP = DIR_SEP_POSIX
        character(:, SK), allocatable   :: errmsg_def
        logical(LK) :: failed_def
        type(shellis_type) :: shellis

        if (present(errmsg)) then
            allocate(errmsg_def, source = errmsg)
        else
            allocate(character(1023, SK) :: errmsg_def)
        end if

        shellis = shellis_type(failed_def, errmsg_def)

        ! Get the Home directory

        failed_def = isFailedGetEnvVar(SK_"HOME", dirHome, errmsg_def) ! Assume Posix shell or PowerShell Core.
        failed_def = failed_def .or. logical(len_trim(dirHome, IK) == 0_IK, LK)
        if (failed_def) then ! assume Windows PowerShell or CMD.
            ! LCOV_EXCL_START
            if (shellis%cmd .or. shellis%powershell) then
                failed_def = isFailedGetEnvVar(SK_"USERPROFILE", dirHome, errmsg_def) ! Windows PowerShell and CMD
                failed_def = failed_def .or. logical(len_trim(dirHome, IK) == 0_IK, LK)
                if (failed_def) then
                    failed_def = isFailedGetEnvVar(SK_"HOMEDRIVE", homeDrive, errmsg_def) ! Windows PowerShell and CMD
                    failed_def = failed_def .or. logical(len_trim(homeDrive, IK) == 0_IK, LK)
                    if (.not. failed_def) then
                        failed_def = isFailedGetEnvVar(SK_"HOMEPATH", homePath, errmsg_def) ! Windows PowerShell and CMD
                        failed_def = failed_def .or. logical(len_trim(homePath, IK) == 0_IK, LK)
                        if (.not. failed_def) dirHome = homeDrive//homePath
                    end if
                end if
            end if
            ! LCOV_EXCL_STOP
        end if

        ! If successful, replace username with the requested user name.

        if (.not. failed_def) then
            failed_def = .not. isDir(dirHome)
            if (failed_def) then
                errmsg_def(:) = SK_"The identified current user home directory does not exist: """//getStr(dirHome)//SK_""""
            elseif (present(user)) then
                !failed_def = logical(len_trim(user) == 0, LK)
                !if (failed_def) then
                !    errmsg_def(:) = SK_"The input user name must not be empty."
                !else
                if (len_trim(user, IK) > 0_IK) then ! substitute username only if the requested username is non-empty.
                    if (shellis%posix .or. shellis%fish) then
                        dirHome = getDirName(dirHome, PSPA)//PSP//user
                    elseif (shellis%windows) then
                        dirHome = getDirName(dirHome, PSWA)//PSW//user
                    end if
                end if
            end if
            !write(*,*) "dirHome"
            !write(*,*) "dirHome"
            !write(*,*) dirHome, failed_def, present(user), isDir(dirHome)
            !if (present(user)) write(*,*) user
        elseif (shellis%posix .or. shellis%fish) then ! try another way on posix shells (including PowerShell Core) and Fish.
            if (present(user)) then ! Fetch the username first.
                user_def = user
            else
                failed_def = isFailedGetEnvVar(SK_"USER", user_def, errmsg_def)
                failed_def = failed_def .or. logical(len_trim(user_def) == 0, LK)
                if (failed_def) then
                    failed_def = isFailedGetEnvVar(SK_"USERNAME", user_def, errmsg_def)
                    failed_def = failed_def .or. logical(len_trim(user_def, IK) == 0_IK, LK)
                    if (failed_def) failed_def = isFailedGetOutput(SK_"whoami", user_def, errmsg_def)
                end if
            end if
            failed_def = failed_def .or. logical(len_trim(user_def, IK) == 0_IK, LK)
            if (.not. failed_def) then
                failed_def = isFailedGetOutput(SK_'getent passwd "'//user_def//SK_'" | cut -d: -f6', dirHome, errmsg_def) ! Get the user-specific path.
                failed_def = failed_def .or. logical(len_trim(dirHome, IK) == 0_IK, LK)
            end if
        end if

        if (failed_def) then
            dirHome = SK_"" ! LCOV_EXCL_LINE
            errmsg_def = MODULE_NAME//SK_"@getDirHome(): Failed to fetch the home directory. "//trim(errmsg_def) ! LCOV_EXCL_LINE
        end if

        if (present(errmsg)) errmsg = errmsg_def ! LCOV_EXCL_LINE

        if (present(failed)) then
            failed = failed_def
        elseif (failed_def) then
            error stop trim(errmsg_def) ! LCOV_EXCL_LINE
        end if

        deallocate(errmsg_def) ! gfortran 11 automatic heap deallocation bug.

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathExpandedUser
        character(2047, SK) :: errmsg
        logical(LK) :: failed
        pathExpandedUser = getPathExpandedUser(path, failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@getPathExpandedUser(): "//trim(errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathExpandedUserFailed
        character(2047, SK) :: errmsg
        pathExpandedUser = getPathExpandedUser(path, failed, errmsg)
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathExpandedUserFailedMsg
        use pm_kind, only: SKC => SK
        use pm_sysShell, only: shellis_type
        character(*,SKC), parameter :: PSPA = DIR_SEP_POSIX_ALL
        character(*,SKC), parameter :: PSWA = DIR_SEP_WINDOWS_ALL
        type(shellis_type) :: shellis
        integer(IK) :: lenPath
        integer(IK) :: firstDirSepPos
        failed = .false._LK
        lenPath = len(path, IK)
        if (lenPath > 0_IK) then
            if (path(1:1) == SKC_"~") then
                shellis = shellis_type(failed, errmsg)
                if (.not. failed) then
                    if (shellis%posix .or. shellis%fish) then
                        firstDirSepPos = scan(path(1:lenPath), PSPA, kind = IK)
                    else
                        firstDirSepPos = scan(path(1:lenPath), PSWA, kind = IK)
                    end if
                    if (firstDirSepPos == 0_IK) firstDirSepPos = lenPath + 1_IK
                    pathExpandedUser = getDirHome(path(2:firstDirSepPos-1), failed = failed, errmsg = errmsg) // path(firstDirSepPos:lenPath)
                    return
                end if
            end if
        end if
        pathExpandedUser = path
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedChangeDir
        ! \warning : The input `dir` must be of default character kind.
#if     INTEL_ENABLED
        use ifport, only: chdir
#endif
#if     INTEL_ENABLED || GNU_ENABLED
        ! This **could** yield an ICE with gfortran on macOS like the other above.
        ! failed = logical(chdir(dir) /= 0, LK)
        failed = chdir(dir) /= 0
#else
        use iso_c_binding, only: c_null_char, c_int
#if     WINDOWS_ENABLED
        interface ! See WIN32 API for details of `_chdir`, `_wchdir`. `chdir` is deprecated
            function chdir(path) result(stat) bind(C, name = "_chdir")
                use iso_c_binding, only: c_char, c_int
                character(1, c_char), intent(in) :: path(*)
                integer(c_int) :: stat
            end function
        end interface
#else
        interface
            function chdir(path) result(stat) bind(C, name = "chdir")
                use iso_c_binding, only: c_char, c_int
                character(1,c_char), intent(in) :: path(*)
                integer(c_int) :: stat
            end function
        end interface
#endif
        failed = logical(chdir(dir//c_null_char) /= int(0,c_int), LK)
#endif
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathTemp
        use pm_sysShell, only: isFailedGetDirTemp
        if (isFailedGetDirTemp(pathTemp) .and. allocated(pathTemp)) deallocate(pathTemp)
        pathTemp = getPathNew(pathTemp, prefix, sep, ext, pid, failed)
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathNew

        use pm_kind, only: SKC => SK
        use pm_val2str, only: getStr
        use pm_sysShell, only: shellis_type
        use pm_parallelism, only: getImageID

        character(1,SKC), parameter     :: USC = SKC_"_"
        character(*,SKC), parameter     :: PSP = DIR_SEP_POSIX
        character(*,SKC), parameter     :: PSW = DIR_SEP_WINDOWS
        character(*,SKC), parameter     :: PSPA = DIR_SEP_POSIX_ALL
        character(*,SKC), parameter     :: PSWA = DIR_SEP_WINDOWS_ALL
        character(:,SKC), allocatable   :: dir_def, prefix_def, ext_def
        integer(IK)     , parameter     :: COUNTER_MAX = 1000_IK
        integer(IK)                     :: counter, lenDir
        integer(IK)                     :: pid_def
        logical(LK)                     :: exists
        type(shellis_type)              :: shellis
        character(8 ,SKC)               :: date
        character(10,SKC)               :: time


        if (present(failed)) then
            shellis = shellis_type(failed)
            if (failed) then
                pathNew = SKC_"" ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            end if
        else
            shellis = shellis_type()
        end if

        if (present(dir)) then
            if (len_trim(dir, IK) > 0_IK) then
                lenDir = len(dir, IK)
                if (shellis%windows .and. scan(PSWA, dir(lenDir:lenDir), kind = IK) == 0_IK) then
                    dir_def = dir//PSW
                elseif (scan(PSPA, dir(lenDir:lenDir), kind = IK) == 0_IK) then
                    dir_def = dir//PSP
                else
                    dir_def = dir
                end if
            else
                dir_def = SKC_""
            end if
        else
            dir_def = SKC_""
        end if

        if (present(prefix)) then
            prefix_def = prefix
        else
            if (present(sep)) then
                prefix_def = SKC_"new"//sep
            else
                prefix_def = SKC_"new"//USC
            end if
        end if

        if (present(ext)) then
            ext_def = ext
        else
            ext_def = SKC_""
        end if

        if (present(pid)) then
            pid_def = pid
        else
            pid_def = getImageID()
        end if


        do counter = 1_IK, COUNTER_MAX
            call date_and_time(date, time)
            if (present(sep)) then
                pathNew = dir_def//prefix_def//date//sep//time(1:6)//sep//time(8:10)//sep//SKC_"pid"//sep//getStr(pid_def)//ext_def
            else
                pathNew = dir_def//prefix_def//date//USC//time(1:6)//USC//time(8:10)//USC//SKC_"pid"//USC//getStr(pid_def)//ext_def
            end if
            inquire(file = pathNew, exist = exists)
            if (.not. exists) return
        end do

        if (present(failed)) then
            failed = .true._LK
            pathNew = SKC_""
        else
            error stop SK_"Failed to generate new file, even after "//getStr(COUNTER_MAX)//SK_" tries."
        end if

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathJoined
        use pm_kind, only: SKC => SK
        character(255, SK):: errmsg
        logical(LK) :: failed
        pathJoined = getPathJoined(p1, p2, failed, errmsg)
        if (failed) error stop MODULE_NAME//SK_"@getPathJoined(): "//trim(errmsg) ! LCOV_EXCL_LINE
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathJoinedFailed

        use pm_kind, only: SKC => SK
        use pm_sysShell, only: shell_type
        character(*,SKC), parameter :: PSWA = DIR_SEP_WINDOWS_ALL
        character(*,SKC), parameter :: PSPA = DIR_SEP_POSIX_ALL
        type(shell_type) :: Shell
        integer(IK) :: lenp1
        integer(IK) :: lenp2

        lenp1 = len(p1, IK)
        lenp2 = len(p2, IK)

        if (lenp1 == 0_IK .and. lenp2 == 0_IK) then
            pathJoined = SKC_""
            failed = .false._LK
            return
        end if

        if (present(errmsg)) then
            shell = shell_type(failed, errmsg)
        else
            shell = shell_type(failed)
        end if
        if (failed) return ! LCOV_EXCL_LINE

        if (shell%is%windows) then
            if (isPathAbsWindows(p2) .or. hasDriveLetter(p2)) then
                pathJoined = p2
            elseif (scan(PSWA, p1(lenp1:lenp1), kind = IK) == 0_IK) then
                pathJoined = p1//shell%dirsep//p2
            else
                pathJoined = p1//p2
            end if
        elseif (shell%is%posix .or. shell%is%fish) then
            if (isPathAbsPosix(p2)) then
                pathJoined = p2
            elseif (scan(PSPA, p1(lenp1:lenp1), kind = IK) == 0_IK) then
                pathJoined = p1//shell%dirsep//p2
            else
                pathJoined = p1//p2
            end if
        end if

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedMakeDir
        use pm_kind, only: SKC => SK

        character(:,SKC), allocatable   :: command
        type(shellis_type)              :: shellis
        integer(IK)                     :: itry
        integer(IK)                     :: ntry_def
        character(*,SKC), parameter     :: LF = new_line(SKC_"a")

        if (present(errmsg)) then
            shellis = shellis_type(failed, errmsg)
        else
            shellis = shellis_type(failed)
        end if

        ! The order of the following conditionals is significant.

        if (shellis%powershell .and. shellis%windows) then
            !>  PowerShell accepts `-p` flag on both POSIX and Windows systems.
            !>  However, in POSIX OS, PowerShell requires `-p` for nested directories.
            command = SKC_"New-Item -ItemType Directory -Path "//getPathVerbatimPowerShell(path)//SKC_" >null 2>&1" ! LCOV_EXCL_LINE
        elseif (shellis%posix) then ! including microsoft powershell core.
            command = SKC_"mkdir -p "//getPathVerbatimPosix(path)//SKC_" >/dev/null 2>&1" ! -p enables nested mkdir
        elseif (shellis%fish) then
            command = SKC_"mkdir -p "//getPathVerbatimFish(path)//SKC_" >/dev/null 2>&1" ! -p enables nested mkdir
        elseif (shellis%cmd) then
            command = SKC_"mkdir "//getPathVerbatimCMD(path)//SKC_" >nul 2>&1" ! LCOV_EXCL_LINE
        else
            !command = SKC_"mkdir "//path ! last shot in the dark. ! LCOV_EXCL_LINE
            if (present(errmsg)) errmsg = MODULE_NAME//SK_"@isFailedMakeDir(): Failed to fetch the shell type."//LF//trim(errmsg) ! LCOV_EXCL_LINE
            failed = .true._LK
            return
        end if
        ! Create the folder.

        if (present(ntry)) then
            CHECK_ASSERTION(__LINE__, ntry > 0_IK, SK_"@isFailedMakeDir(): The input `ntry` must be a positive number. ntry = "//getStr(ntry)) ! fpp
            ntry_def = ntry
        else
            ntry_def = 1_IK
        end if

        do itry = 1_IK, ntry_def
            failed = isFailedExec(command, wait = wait, cmdmsg = errmsg)
            failed = .not. isDir(path)
            if (.not. failed) return
        end do

        if (failed .and. present(errmsg)) errmsg = MODULE_NAME//SK_"@isFailedMakeDir(): Failed to create the requested directory."//LF//trim(errmsg) ! LCOV_EXCL_LINE
        ! \warning do NOT remove. \bug gfortran with heap memory allocations fails to automatically deallocate on exit.
        deallocate(command)

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure isFailedMakeDirTemp

        use pm_kind, only: SKC => SK
        use pm_sysShell, only: isFailedGetDirTemp
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedMakeDirTemp()"
        integer(IK)     , parameter :: NTRY = 3_IK
        character(:,SKC), allocatable :: parent_def
        integer(IK) :: itry

        do itry = 1_IK, NTRY

            if (present(parent)) then
                path = getPathNew(dir = parent, failed = failed)
            elseif (.not. isFailedGetDirTemp(parent_def)) then
                path = getPathNew(dir = parent_def, prefix = SK_"temp_", sep = SK_"_", failed = failed)
            else
                path = getPathNew(prefix = SK_"temp_", sep = SK_"_", failed = failed)
            end if

            if (failed) then
                if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_"Failed to generate a unique dir name via getPathNew()." ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            end if

            if (.not. isFailedMakeDir(path, wait = .true._LK, ntry = NTRY, errmsg = errmsg)) return

        end do

        failed = .true._LK
        if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg) ! LCOV_EXCL_LINE

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedCopy_ENABLED 1

    module procedure isFailedCopy
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef isFailedCopy_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedMove_ENABLED 1

    module procedure isFailedMove
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef isFailedMove_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedRemove_ENABLED 1

    module procedure isFailedRemove
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

#undef isFailedRemove_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure glob_BSSK
        use pm_kind, only: SKC => SK
        character(511, SK) :: errmsg
        if (isFailedGlob(pattern, list, errmsg = errmsg)) error stop MODULE_NAME//SK_"@glob(): "//trim(errmsg)
    end procedure

#define isFailedGlob_ENABLED 1

#define SK_ENABLED 1
    module procedure isFailedGlob_SK
        use pm_kind, only: SKC => SK
        use pm_strASCII, only: setAsciiFromEscaped
        use pm_sysShell, only: isFailedGetOutput, shell_type
#include "pm_sysPath@routines.inc.F90"
    end procedure
#undef SK_ENABLED

#define BSSK_ENABLED 1
    module procedure isFailedGlob_BSSK
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure
#undef BSSK_ENABLED

#undef isFailedGlob_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure ls_BSSK
        use pm_kind, only: SKC => SK
        character(511, SK) :: errmsg
        if (isFailedList(path, list, sort = sort, showdir = showdir, showfile = showfile, showhidden = showhidden, reversed = reversed, errmsg = errmsg)) error stop MODULE_NAME//SK_"@ls(): "//trim(errmsg)
    end procedure

#define isFailedList_ENABLED 1

#define SK_ENABLED 1
    module procedure isFailedList_SK
        use pm_kind, only: SKC => SK
        use pm_strASCII, only: setAsciiFromEscaped
        use pm_sysShell, only: isFailedGetOutput, shell_type
#include "pm_sysPath@routines.inc.F90"
    end procedure
#undef SK_ENABLED

#define BSSK_ENABLED 1
    module procedure isFailedList_BSSK
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure
#undef BSSK_ENABLED

#undef isFailedList_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define hasDriveLetter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure hasDriveLetter_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure hasDriveLetter_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure hasDriveLetter_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure hasDriveLetter_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure hasDriveLetter_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef hasDriveLetter_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isPathAbsPosix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure isPathAbsPosix_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isPathAbsPosix_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isPathAbsPosix_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isPathAbsPosix_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isPathAbsPosix_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isPathAbsPosix_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isPathAbsWindows_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure isPathAbsWindows_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isPathAbsWindows_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isPathAbsWindows_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isPathAbsWindows_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isPathAbsWindows_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isPathAbsWindows_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDirSep_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirSep_SK
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getDirSep_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getDirSep_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getDirSep_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getDirSep_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getDirSep_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDirSep_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDirSeps_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getDirSeps_SK
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getDirSeps_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getDirSeps_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getDirSeps_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getDirSeps_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getDirSeps_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDirSeps_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathSep_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getPathSep_SK
        use pm_kind, only: SKC => SK
#include "pm_sysPath@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getPathSep_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathSep_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathSep_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathSep_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathSep_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathSep_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathVerbatimCMD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getPathVerbatimCMD_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathVerbatimCMD_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathVerbatimCMD_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathVerbatimCMD_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathVerbatimCMD_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathVerbatimCMD_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathVerbatimPowerShell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getPathVerbatimPowerShell_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathVerbatimPowerShell_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathVerbatimPowerShell_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathVerbatimPowerShell_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathVerbatimPowerShell_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathVerbatimPowerShell_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathVerbatimPosix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getPathVerbatimPosix_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathVerbatimPosix_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathVerbatimPosix_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathVerbatimPosix_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathVerbatimPosix_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathVerbatimPosix_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathVerbatimFish_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getPathVerbatimFish_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathVerbatimFish_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathVerbatimFish_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathVerbatimFish_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathVerbatimFish_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathVerbatimFish_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathPosix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathPosix_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPathPosix_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathPosix_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathPosix_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathPosix_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathPosix_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#undef getPathPosix_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathPosix_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathPosix_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathPosix_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPathPosix_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPathPosix_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPathPosix_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPathPosix_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPathPosix_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#undef setPathPosix_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPathPosix_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathPosixEscaped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathPosixEscaped_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPathPosixEscaped_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathPosixEscaped_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathPosixEscaped_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathPosixEscaped_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathPosixEscaped_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#undef getPathPosixEscaped_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathPosixEscaped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathPosixEscaped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathPosixEscaped_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPathPosixEscaped_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPathPosixEscaped_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPathPosixEscaped_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPathPosixEscaped_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPathPosixEscaped_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#undef setPathPosixEscaped_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPathPosixEscaped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathWindows_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathWindows_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure getPathWindows_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathWindows_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathWindows_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathWindows_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathWindows_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#undef getPathWindows_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathWindows_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathWindows_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPathWindows_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure setPathWindows_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setPathWindows_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setPathWindows_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setPathWindows_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setPathWindows_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#undef setPathWindows_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPathWindows_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPathHostNameIndex_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module procedure getPathHostNameIndex_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getPathHostNameIndex_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getPathHostNameIndex_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getPathHostNameIndex_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getPathHostNameIndex_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sysPath@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPathHostNameIndex_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getPathNew_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getPathNew_D0_SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getPathNew_D0_SK5
!        use pm_kind, only: SKC => SK5
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getPathNew_D0_SK4
!        use pm_kind, only: SKC => SK4
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getPathNew_D0_SK3
!        use pm_kind, only: SKC => SK3
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getPathNew_D0_SK2
!        use pm_kind, only: SKC => SK2
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getPathNew_D0_SK1
!        use pm_kind, only: SKC => SK1
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#undef getPathNew_D0_SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getPathNew_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define isDir_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define isDir_D0_SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure isDir_D0_SK5
!        use pm_kind, only: SKC => SK5
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure isDir_D0_SK4
!        use pm_kind, only: SKC => SK4
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure isDir_D0_SK3
!        use pm_kind, only: SKC => SK3
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure isDir_D0_SK2
!        use pm_kind, only: SKC => SK2
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure isDir_D0_SK1
!        use pm_kind, only: SKC => SK1
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#undef isDir_D0_SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef isDir_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define isFile_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define isFile_D0_SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure isFile_D0_SK5
!        use pm_kind, only: SKC => SK5
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure isFile_D0_SK4
!        use pm_kind, only: SKC => SK4
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure isFile_D0_SK3
!        use pm_kind, only: SKC => SK3
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure isFile_D0_SK2
!        use pm_kind, only: SKC => SK2
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure isFile_D0_SK1
!        use pm_kind, only: SKC => SK1
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#undef isFile_D0_SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef isFile_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define mkdir_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define mkdir_D0_SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure mkdir_D0_SK5
!        use pm_kind, only: SKC => SK5
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure mkdir_D0_SK4
!        use pm_kind, only: SKC => SK4
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure mkdir_D0_SK3
!        use pm_kind, only: SKC => SK3
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure mkdir_D0_SK2
!        use pm_kind, only: SKC => SK2
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure mkdir_D0_SK1
!        use pm_kind, only: SKC => SK1
!#include "pm_sysPath@routines.inc.F90"
!    end procedure
!#endif
!
!#undef mkdir_D0_SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef mkdir_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines

#undef  UNESCAPE_REINDEX_LIST
#undef  CHECK_ASSERTION
#undef  ERROR_STOP_IF
