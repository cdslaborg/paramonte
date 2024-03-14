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
!>  This module contains implementations of the procedures in [pm_sysPath](@ref pm_sysPath).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     INTEL_ENABLED && WINDOWS_ENABLED
#define SHARED, shared
#else
#define SHARED
#endif
        !%%%%%%%%%%%%
#if     isDir_ENABLED
        !%%%%%%%%%%%%

        character(*,SKC), parameter :: DSP = DIR_SEP_POSIX
        character(*,SKC), parameter :: DSW = DIR_SEP_WINDOWS
#if     INTEL_ENABLED && DD_ENABLED
        inquire(directory = path, exist = pathIsDir)
#elif   INTEL_ENABLED && II_ENABLED
        if (present(iomsg)) then
            inquire(directory = path, exist = pathIsDir, iostat = iostat, iomsg = iomsg)
        else
            inquire(directory = path, exist = pathIsDir, iostat = iostat)
        end if
#elif   GNU_ENABLED && DD_ENABLED
        inquire(file = path//DSP, exist = pathIsDir)
        if (.not. pathIsDir) then
            inquire(file = path//DSW, exist = pathIsDir)
        end if
#elif   GNU_ENABLED && II_ENABLED
        if (present(iomsg)) then
            inquire(file = path//DSP, exist = pathIsDir, iostat = iostat, iomsg = iomsg)
        else
            inquire(file = path//DSP, exist = pathIsDir, iostat = iostat)
        end if
#elif   DD_ENABLED || II_ENABLED
        character(len(c_null_char),SKC), parameter :: NULL_CHAR = c_null_char
        character(len(path) + 1, c_char) :: pathc
        interface
            function pm_sys_isdirc(pathc) result(itis) bind(C, name = "pm_sys_isdirc")
                import :: c_int32_t, c_char
                character(1, c_char), intent(in) :: pathc(*)
                integer(c_int32_t) :: itis
            end function
        end interface
        pathc = path//NULL_CHAR
        pathIsDir = logical(pm_sys_isdirc(pathc) == 1_c_int32_t, LK)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%
#elif   isFile_ENABLED
        !%%%%%%%%%%%%%

#if     DD_ENABLED
        inquire(file = path, exist = pathIsFile)
#elif   II_ENABLED
        if (present(iomsg)) then
            inquire(file = path, exist = pathIsFile, iostat = iostat, iomsg = iomsg)
        else
            inquire(file = path, exist = pathIsFile, iostat = iostat)
        end if
        if (iostat /= 0_IK) return
#else
#error  "Unrecognized interface."
#endif
#if     !INTEL_ENABLED
        if (pathIsFile) then
#if         DD_ENABLED
            pathIsFile = .not. isDir(path)
#elif       II_ENABLED
            pathIsFile = .not. isDir(path, iostat, iomsg)
#endif
        end if
#endif

        !%%%%%%%%%%%%%%%
#elif   isExtant_ENABLED
        !%%%%%%%%%%%%%%%

#if     DD_ENABLED
        extant = isFile(path)
        if (.not. extant) extant = isDir(path)
#elif   II_ENABLED
        extant = isFile(path, iostat, iomsg)
        if (.not. extant) extant = isDir(path, iostat, iomsg)
#else
#error  "Unrecognized interface."
#endif
!#if     DD_ENABLED
!        inquire(file = path, exist = extant)
!#elif   II_ENABLED
!        if (present(iomsg)) then
!            inquire(file = path, exist = extant, iostat = iostat, iomsg = iomsg)
!        else
!            inquire(file = path, exist = extant, iostat = iostat)
!        end if
!        if (iostat /= 0_IK) return
!#else
!#error  "Unrecognized interface."
!#endif
!#if     !INTEL_ENABLED
!        if (.not. extant) then
!#if         DD_ENABLED
!            extant = isDir(path)
!#elif       II_ENABLED
!            extant = isDir(path, iostat, iomsg)
!#endif
!        end if
!#endif

        !%%%%%%%%%%%%%%%%%
#elif   getDirName_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: index
#if     Def_ENABLED
        index = getIndexDirName(path, dirsep)
        if (0_IK < index) then
            dirname = path(1_IK : index)
        else
            dirname = SKC_"."
        end if
#elif   PM_ENABLED
        index = getIndexDirName(path, dirsep, style)
        dirname = path(1_IK : index)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%
#elif   getBaseName_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     Def_ENABLED
        if (1_IK < len(path, IK)) then
            block
                integer(IK) :: offset
                ! remove trailing directory separator duplicates.
                offset = len(path, IK)
                if (len(dirsep, IK) == 1_IK) then
                    do
                        if (offset < 2_IK) exit
                        if (dirsep /= path(offset : offset)) exit
                        offset = offset - 1_IK
                    end do
                else
                    do
                        if (offset < 2_IK) exit
                        if (index(dirsep, path(offset : offset), kind = IK) == 0_IK) exit
                        offset = offset - 1_IK
                    end do
                end if
                basename = path(getIndexBaseName(path(1 : offset), dirsep) : offset)
            end block
        else
            basename = path
        end if
#elif   PM_ENABLED
        basename = path(getIndexBaseName(path, dirsep, style) :)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getIndexDirName_ENABLED || getIndexBaseName_ENABLED) && PM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0_IK < len(dirsep, IK), MODULE_NAME//SK_": The input `0 < len(dirsep)` must be a positive number. len(dirsep) = "//getStr(len(dirsep, IK)))
        index = scan(path, dirsep, back = .true., kind = IK)
#if     getIndexBaseName_ENABLED
        index = index + 1_IK
#elif   !getIndexDirName_ENABLED
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getIndexDirName_ENABLED || getIndexBaseName_ENABLED) && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: offset
        CHECK_ASSERTION(__LINE__, 0_IK < len(dirsep, IK), MODULE_NAME//SK_": The input `0 < len(dirsep)` must be a positive number. len(dirsep) = "//getStr(len(dirsep, IK)))
#if     getIndexDirName_ENABLED
        if (len(path, IK) < 2_IK) then
            index = scan(path, dirsep, back = .true., kind = IK)
            !index = index + 1_IK
            return
        end if
#endif
        ! remove trailing directory separator duplicates.
        offset = len(path, IK)
        if (len(dirsep, IK) == 1_IK) then
            do
                if (offset < 2_IK) exit
                if (dirsep /= path(offset : offset)) exit
                offset = offset - 1_IK
            end do
        else
            block
                intrinsic :: index
                do
                    if (offset < 2_IK) exit
                    if (index(dirsep, path(offset : offset), kind = IK) == 0_IK) exit
                    offset = offset - 1_IK
                end do
            end block
        end if
#if     getIndexDirName_ENABLED
        block
            intrinsic :: index
            offset = scan(path(1 : offset), dirsep, back = .true., kind = IK)
            do
                if (offset < 2_IK) exit
                offset = offset - 1_IK
                if (index(dirsep, path(offset : offset), kind = IK) == 0_IK) exit
            end do
            if (offset == 2_IK) then ! `$(dirname "../")` on POSIX command line yields `"."`
                if (path(1:2) == SKC_"..") offset = 1_IK
            end if
        end block
        index = offset
#elif   getIndexBaseName_ENABLED
        index = scan(path(1 : offset - 1), dirsep, back = .true., kind = IK) + 1_IK
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDirSep_ENABLED || getDirSeps_ENABLED || getPathSep_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! \warning
        ! We cannot use FPP macros because this is dependent on the runtime shell.
        type(shell_type) :: shell
        if (present(failed) .and. present(errmsg)) then
            shell = shell_type(failed, errmsg)
        elseif (present(failed)) then
            shell = shell_type(failed)
        else
            shell = shell_type()
        end if
#if     getDirSep_ENABLED
        dirsep = shell%dirsep
#elif   getDirSeps_ENABLED
        dirseps = shell%dirseps
#elif   getPathSep_ENABLED
        pathsep = shell%pathsep
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPathVerbatimCMD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        pathVerbatim = SKC_'"'//getRemoved(path, SKC_'"')//SKC_'"'

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPathVerbatimPowerShell_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pathVerbatim = SKC_"'"//getReplaced(path, SKC_"'", SKC_"''")//SKC_"'"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPathVerbatimPosix_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pathVerbatim = SKC_"'"//getReplaced(path, SKC_"'", SKC_"'\''")//SKC_"'"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPathVerbatimFish_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pathVerbatim = SKC_"'"//getReplaced(getReplaced(path, SKC_"\", SKC_"\\"), SKC_"'", SKC_"\'")//SKC_"'"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPathPosixEscaped_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setPathPosixEscaped(pathEscaped, path)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPathPosixEscaped_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(1,SKC), parameter :: POSIX_RESERVED_CHR_SKC(*) = POSIX_RESERVED_CHR
        integer(IK) :: i, j, counter, lenPath, Loc(len(path, IK))
        lenPath = len(path, IK)
        if (lenPath == 0_IK) then
            pathEscaped = path
        else
            counter = 0_IK
            do i = 1_IK, lenPath
                loopOverShellEscapeChars: do j = 1_IK, size(POSIX_RESERVED_CHR_SKC, 1, IK)
                    if (path(i:i) == POSIX_RESERVED_CHR_SKC(j)) then
                        counter = counter + 1_IK
                        Loc(counter) = i
                        exit loopOverShellEscapeChars
                    end if
                end do loopOverShellEscapeChars
            end do
            allocate(character(lenPath + counter, SKC) :: pathEscaped)
            call setInserted(pathEscaped, path, SKC_"\", Loc(1:counter))
        end if

        !%%%%%%%%%%%%%%%%%%%
#elif   getPathPosix_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        pathPosix = path
        call setPathPosix(pathPosix, ignore)

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getPathWindows_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        pathWindows = path
        call setPathWindows(pathWindows, ignore)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPathPosix_ENABLED || setPathWindows_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(:,SKC), allocatable   :: pathNew
        character(*,SKC), parameter     :: PSA = DIR_SEP_ALL ! Path Separators All
#if     setPathPosix_ENABLED
        character(1,SKC), parameter     :: PSN = DIR_SEP_POSIX ! Path Separator Native
#elif   setPathWindows_ENABLED
        character(1,SKC), parameter     :: PSN = DIR_SEP_WINDOWS ! Path Separator Native
#endif
        integer(IK)                     :: i, j, lenPath, lenPathMin, lenIgnore, lenPathNew!, lenSeg, istart
       !logical(LK)                     :: isDot ! is pattern `PSN//"."//PSN`.
        logical(LK)                     :: isSPS ! is significant slash (meaning that it is not preceded by a slash).

        !   \warning
        !   The following code assumes that the Windows and posix directory separators will not change for the foreseeable future.
        !   As a layer of safety, the following code checks for the length DSWA to be 2 and its contents to have the posix directory separator.
        !error_stop_if(len(DSWA) /= 2, SK_"@setPathPosix()/setPathWindows(): Internal library error occurred. The condition `len(DSWA) == 2` must hold.")
        !error_stop_if(len(DSWA) /= 2, SK_"@setPathPosix()/setPathWindows(): Internal library error occurred. The condition `len(DSWA) == 2` must hold.")
        !error_stop_if(index(DSWA,PSN) == 0, SK_"@setPathPosix()/setPathWindows(): Internal library error occurred. The condition `index(DSWA,PSN) > 0` must hold.")

        if (present(ignore)) then
            lenIgnore = len(ignore, IK)
        else
            lenIgnore = 0_IK
        end if

        lenPath = len(path, IK)
        if (lenPath == 0_IK) return

        !   Note that multiple `\` characters in sequence in Linux or Windows reduce to a single `\`. Therefore, ignore all `\` characters.
        !   Here we are assuming that the input path is convertible to a valid Windows path, meaning that no Windows reserved character appears in the path.
        !   If there are any such illegal characters (for example, double quotation marks) or Windows reserved words, we do not check for them here.
        !   Illegal Windows characters (particularly the double quotation mark) can be handled via `getPathVerbatimCMD()` and `getPathVerbatimPowerShell()`.

        allocate(character(lenPath,SKC) :: pathNew)

        ! Skip any UNC path root.

        lenPathNew = getPathHostNameIndex(path, PSA)
        if (lenPathNew > 0_IK) then
            if (lenIgnore == 0_IK) then
                pathNew(1:2) = PSN//PSN
            elseif (lenIgnore == 2_IK .and. pathNew(1:2) == ignore) then
                pathNew(1:2) = ignore
            else
                pathNew(1:2) = PSN//PSN
            end if
            pathNew(3:lenPathNew) = path(3:lenPathNew)
        end if

        ! Construct the normalized path.

        isSPS = .true._LK
       !isSPS = .false._LK
       !isDot = .false._LK
       !istart = lenPathNew + 1_IK
       !i = istart
        i = lenPathNew + 1_IK
        lenPathMin = max(1_IK, lenPathNew)
        if (lenIgnore == 0_IK) then
            ! Remove trailing directory separators.
            !do
            !    if (lenPath == lenPathMin) exit
            !    if (scan(PSA, path(lenPath:lenPath), kind = IK) == 0_IK) exit
            !    lenPath = lenPath - 1_IK
            !end do
            do
                if (i > lenPath) exit
!#define         INCREMENT_PATH_NEW \
!                if (scan(PSA, path(i:i), kind = IK) == 0_IK) then; \
!                    lenPathNew = lenPathNew + 1_IK; \
!                    !isDot = path(i:i) == SKC_"." .and. .not. isSPS; \
!                    pathNew(lenPathNew:lenPathNew) = path(i:i); \
!                    isSPS = .true._LK; \
!                elseif (isSPS) then; \
!               !elseif (isSPS .or. i == istart) then; \
!                    !if (isDot .and. lenPath > 2_IK) then; \
!                    !    lenPathNew = lenPathNew - 1_IK; \
!                    !else; \
!                        lenPathNew = lenPathNew + 1_IK; \
!                        pathNew(lenPathNew:lenPathNew) = PSN; \
!                    !end if; \
!                    isSPS = .false._LK; \
!                end if; \
!                i = i + 1_IK;
#define         INCREMENT_PATH_NEW \
                if (scan(PSA, path(i:i), kind = IK) == 0_IK) then; \
                    lenPathNew = lenPathNew + 1_IK; \
                    pathNew(lenPathNew:lenPathNew) = path(i:i); \
                    isSPS = .true._LK; \
                elseif (isSPS) then; \
                    lenPathNew = lenPathNew + 1_IK; \
                    pathNew(lenPathNew:lenPathNew) = PSN; \
                    isSPS = .false._LK; \
                end if; \
                i = i + 1_IK;
                INCREMENT_PATH_NEW
            end do
        else
            !! first start index.
            !j = (lenPath / lenIgnore) * lenIgnore
            !if (j == lenPath) then
            !    j = lenPath - lenIgnore + 1_IK
            !else
            !    j = min((lenPath / lenIgnore) * lenIgnore + 1_IK, lenPath)
            !end if
            !! Trim dirsep from the end.
            !do
            !    if (lenPath == lenPathMin) exit
            !    lenSeg = min(lenPath,j+lenIgnore-1_IK) - j + 1_IK
            !    if (lenSeg == lenIgnore .and. path(j:j+lenSeg-1_IK) == ignore) exit
            !    if (scan(PSA, path(lenPath:lenPath), kind = IK) == 0_IK) exit
            !    lenPath = lenPath - 1_IK
            !    j = max(1_IK, j - 1_IK)
            !end do
            do
                if (i > lenPath) exit
                j = min(i + lenIgnore - 1_IK, lenPath)
                if (j - i + 1_IK == lenIgnore .and. path(i:j) == ignore) then
                    lenPathNew = lenPathNew + 1_IK ! Do not move this.
                    pathNew(lenPathNew : lenPathNew + lenIgnore - 1_IK) = ignore
                    lenPathNew = lenPathNew + lenIgnore - 1_IK
                    i = i + lenIgnore
                else
                    INCREMENT_PATH_NEW
                !else ! This should happen only occasionally when `j == lenPath` and path(i:j) must be padded for comparison with `ignore`.
                !    lenPathNew = lenPathNew + 1_IK ! Do not move this.
                !    pathNew(lenPathNew : lenPathNew + j - i + 1_IK) = path(i:j)
                !    lenPathNew = lenPathNew + j - i + 1_IK
                !    i = i + j - i + 2_IK
                end if
            end do
        end if
        !#if     setPathPosix_ENABLED
        !        if (present(escaped)) then
        !            if (escaped) call setPathPosixEscaped(path, pathNew(1:lenPathNew))
        !            deallocate(pathNew) ! \bug : required because gfortran automatic deallocation of heap memory fails as of version 12.
        !            return
        !        end if
        !#endif
        if (lenPathNew > 1_IK .and. pathNew(lenPathNew:lenPathNew) == PSN) lenPathNew = lenPathNew - 1_IK
        path = pathNew(1:lenPathNew)
        deallocate(pathNew) ! \bug : required because gfortran automatic deallocation of heap memory fails as of version 12.

#undef  INCREMENT_PATH_NEW

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPathHostNameIndex_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenPath
        index = 0_IK
        lenPath = len(path, IK)
        if (lenPath < 3_IK) return
        if (scan(dirsep, path(1:1), kind = IK) == 0_IK) return
        if (scan(dirsep, path(2:2), kind = IK) == 0_IK) return
        if (scan(dirsep, path(3:3), kind = IK)  > 0_IK) return
        index = 3_IK
        do
            index = index + 1_IK
            if (index > lenPath) exit
            if (scan(dirsep, path(index:index), kind = IK) > 0_IK) exit
        end do
        index = index - 1_IK

        !%%%%%%%%%%%%%%%%%%%%%
#elif   hasDriveLetter_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        if (len(path,IK) > 1_IK) then
            pathHasDriveLetter = isCharAlpha(path(1:1)) .and. path(2:2) == SKC_":"
        else
            pathHasDriveLetter = .false._LK
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   isPathAbsWindows_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        character(*,SKC), parameter :: DSWA = DIR_SEP_WINDOWS_ALL
        ERROR_STOP_IF(len(DSWA) /= 2, MODULE_NAME//SK_"@isPathAbsWindows(): Internal library error occurred. The condition `len(DSWA) == 2` must hold.")

        if (len(path, IK) > 2_IK) then
            !pathIsAbs = logical(index(DSWA, path(1:1), kind = IK) > 0_IK, LK)
            pathIsAbs = hasDriveLetter(path) .and. scan(DSWA, path(3:3), kind = IK) > 0_IK ! (path(3:3) == DSWA(1:1) .or. path(3:3) == DSWA(2:2)) ! \bug gfortran 11 bug: Error: Operands of comparison operator '==' at (1) are CHARACTER(*,4)/CHARACTER(1)
            if (pathIsAbs) return
        end if

        if (len(path, IK) > 1_IK) then
           !pathIsAbs = (path(1:1) == DSWA(1:1) .or. path(1:1) == DSWA(2:2)) .and. (path(2:2) == DSWA(1:1) .or. path(2:2) == DSWA(2:2)) ! \bug gfortran 11 bug: Error: Operands of comparison operator '==' at (1) are CHARACTER(*,4)/CHARACTER(1)
            pathIsAbs = verify(path(1:2), DSWA, kind = IK) == 0_IK
            if (pathIsAbs) return
        end if

        pathIsAbs = .false._LK

        !%%%%%%%%%%%%%%%%%%%%%
#elif   isPathAbsPosix_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        character(*,SKC), parameter :: DSPA = DIR_SEP_POSIX_ALL
        ERROR_STOP_IF(len(DSPA) > 1, MODULE_NAME//SK_"@isPathAbsPosix(): Internal error occurred. `len(DSPA) == 1` must hold.")
        if (len(path, IK) > 0_IK) then
            pathIsAbs = logical(path(1:1) == DSPA, LK)
        else
            pathIsAbs = .false._LK
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedGlob_ENABLED && BSSK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iloc, nloc
        integer(IK), allocatable :: index(:,:)
        character(:,SKC), allocatable :: contents
        failed = isFailedGlob(pattern, contents, index, errmsg)
        nloc = size(index, 2, IK)
        allocate(list(nloc))
        do concurrent(iloc = 1 : nloc)
            list(iloc)%val = contents(index(1, iloc) : index(2, iloc))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedGlob_ENABLED && SK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: NLC = new_line(SKC_"a")
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedGlob()"
        character(:,SKC), allocatable :: command, errmsg_def
        integer(IK) :: exitstat, iell
        type(shell_type) :: shell

        if (present(errmsg)) then
            allocate(character(len(errmsg, IK), SK) :: errmsg_def)
        else
            allocate(character(31, SK) :: errmsg_def)
        end if

        ! Infer the runtime shell.

        shell = shell_type(failed, errmsg_def)
        if (failed) then
            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg_def) ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        ! Determine sorting method.

        if (shell%is%powershell .or. shell%is%cmd) then
            command = SKC_"(Resolve-Path "//getPathVerbatimPowerShell(pattern)//SKC_").path"
            if (shell%is%cmd) command = SKC_"powershell -command "//command
        elseif (shell%is%bash .or. shell%is%ksh .or. shell%is%zsh) then
            command = getCommandBashStyle(pattern)
        else
            ! check if bash exists on the system.
            failed = isFailedExec(SKC_"bash --version > /dev/null 2>&1")
            if (failed) then
                if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Unsupported runtime shell: "//shell%name ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            end if
            command = SKC_"bash -c "//getPathVerbatimPosix(getCommandBashStyle(pattern))
        end if

        ! Fetch the list.

        failed = isFailedGetOutput(command, list, errmsg_def, exitstat = exitstat)
        if (failed) then
            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg_def) ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        elseif (exitstat /= 0_IK) then
            index = reshape([integer(IK) ::], [0,0])
            list = SKC_""
            return
        end if

        exitstat = 0_IK
        call setSplit(index, list, NLC)
        do iell = 1, size(index, 2, IK)
            if (.not. isExtant(list(index(1, iell) : index(2, iell)), exitstat, errmsg)) then
                index = reshape([integer(IK) ::], [0, 0])
                list = SKC_""
                return
            elseif (exitstat /= 0) then
                failed = .true._LK ! LCOV_EXCL_LINE
                if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Failed while verifying the existence of """// & ! LCOV_EXCL_LINE
                list(index(1, iell) : index(2, iell))//SK_""". "//trim(errmsg) ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            end if
        end do

    contains

        PURE function getCommandBashStyle(pattern) result(command)
            character(*,SKC), intent(in) :: pattern
            character(:,SKC), allocatable :: command
            command = getReplaced(getReplaced(getPathVerbatimPosix(pattern), SKC_"*", SKC_"'*'"), SKC_"?", SKC_"'?'")
            command = SKC_"list=("//command//SKC_"); for file in ""${list[@]}""; do echo ""$file""; done"
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedList_ENABLED && BSSK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iloc, nloc
        integer(IK), allocatable :: index(:,:)
        character(:,SKC), allocatable :: contents
        failed = isFailedList(path, contents, index, sort, showdir, showfile, showhidden, reversed, errmsg)
        nloc = size(index, 2, IK)
        allocate(list(nloc))
        do concurrent(iloc = 1 : nloc)
            list(iloc)%val = contents(index(1, iloc) : index(2, iloc))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedList_ENABLED && SK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: NLC = new_line(SKC_"a")
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedList()"
        character(:,SKC), allocatable :: command, dirlist, filelist, errmsg_def
        integer(IK) , allocatable :: dirIndex(:,:), fileIndex(:,:)
        integer(IK) :: lenList, istart, endloc, i, exitstat
        logical(LK) :: showhidden_def
        logical(LK) :: reversed_def
        logical(LK) :: showfile_def
        logical(LK) :: showdir_def
        type(shell_type) :: shell

        if (present(errmsg)) then
            allocate(character(len(errmsg,IK), SK) :: errmsg_def)
        else
            allocate(character(31, SK) :: errmsg_def)
        end if

        if (present(reversed)) then
            reversed_def = reversed
        else
            reversed_def = .false._LK
        end if

        if (present(showdir)) then
            showdir_def = showdir
        else
            showdir_def = .true._LK
        end if

        if (present(showfile)) then
            showfile_def = showfile
        else
            showfile_def = .true._LK
        end if

        if (present(showhidden)) then
            showhidden_def = showhidden
        else
            showhidden_def = .true._LK
        end if

        ! Infer the runtime shell.

        shell = shell_type(failed, errmsg_def)
        if (failed) then
            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg_def) ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE
        end if

        ! Determine sorting method.

        if (shell%is%posix .or. shell%is%fish) then

            !   `ls` flags:
            !   -a, --all                   :   do not ignore entries starting with `.`.
            !   -Q, --quote-name            :   enclose entry names in double quotes.
            !   -A, --almost-all            :   do not list implied `.` and `..`.
            !   -b, --escape                :   print C-style escapes for nongraphic characters.
            !   -c with -lt                 :   sort by, and show, ctime (time of last modification of file status information); with -l: show ctime and sort by name; otherwise: sort by ctime, newest first.
            !   -F, --classify              :   append indicator (one of */=>@|) to entries. Display
            !                                       -#  a slash (/) immediately after each pathname that is a directory,
            !                                       -#  an asterisk (`*') after each that is executable,
            !                                       -#  an at sign (`@') after each symbolic link,
            !                                       -#  a percent sign (`%') after each whiteout,
            !                                       -#  an equal sign (`=') after each socket,
            !                                       -#  avertical bar (`|') after each that is a FIFO.
            !   --group-directories-first   :   group directories before files; can be augmented with a `--sort` option, but any use of `--sort=none` (`-U`) disables grouping.
            !   -G, --no-group              :   in a long listing, don't print group names
            !   -h, --human-readable        :   with -l and -s, print sizes like 1K 234M 2G etc.
            !   -H                          :   Evaluate the file information and file type for symbolic links specified on the command line to be those of the file referenced by the link, and not the link itself;
            !                                   however, ls shall write the name of the link itself and not the file referenced by the link.
            !   -L                          :   Evaluate the file information and file type for all symbolic links (whether named on the command line or encountered in a file hierarchy)
            !                                   to be those of the file referenced by the link, and not the link itself; however, ls shall write the name of the link itself and not the file referenced by the link.
            !                                   When -L is used with -l, write the contents of symbolic links in the long format (see the STDOUT section).
            !   -k, --kibibytes             :   default to 1024-byte blocks for disk usage; used only with `-s` and per directory totals.
            !   -m                          :   fill width with a comma separated list of entries.
            !   -p, --indicator-style=slash :   append / indicator to directories.
            !   -Q, --quote-name            :   enclose entry names in double quotes.
            !   --quoting-style=WORD        :   use quoting style WORD for entry names: `literal`, `locale`, `shell`, `shell-always`, `shell-escape`, `shell-escape-always`, `c`, `escape` (overrides QUOTING_STYLE environment variable).
            !   -r, --reverse               :   reverse order while sorting.
            !   -S                          :   sort by file size, largest first.
            !   --sort=WORD                 :   sort by WORD instead of name: none (-U), size (-S), time (-t), version (-v), extension (-X).
            !   -X                          :   sort alphabetically by entry extension.
            !   -t                          :   sort by time, newest first.
            !   -1                          :   list one file per line.  Avoid '\n' with -q or -b.
            !   -u                          :   Use time of last access instead of last modification of the file for sorting (-t) or writing (-l).
            !
            !   \see
            !   https://pubs.opengroup.org/onlinepubs/9699919799/utilities/ls.html
            !
            !   \warning
            !   Grouping directories first in the command below is crucial for the proper functioning
            !   of the rest of the code below (when either `showdir` or `showfile` is `.false.`).
            command = SKC_"ls --group-directories-first -bpG"
            if (showhidden_def) command = command//SKC_"A"
            if (present(sort)) then
                if (sort /= SK_"name") then
                    if (reversed_def) command = command//SKC_"r"
                elseif (sort == SK_"tmod") then
                    if (reversed_def) then
                        command = command//SKC_"t"
                    else
                        command = command//SKC_"tr"
                    end if
                elseif (sort == SK_"tacc") then
                    if (reversed_def) then
                        command = command//SKC_"u"
                    else
                        command = command//SKC_"ur"
                    end if
                elseif (sort == SK_"size") then
                    if (reversed_def) then
                        command = command//SKC_"S"
                    else
                        command = command//SKC_"Sr"
                    end if
                elseif (sort == SK_"fext") then
                    if (reversed_def) then
                        command = command//SKC_"Xr"
                    else
                        command = command//SKC_"X"
                    end if
                else
                    if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Unrecognized value for the `sort` input argument: "//sort
                    failed = .true._LK
                    return
                end if
            else
                if (reversed_def) command = command//SKC_" -r"
            end if

            ! Fetch the list.
            failed = isFailedGetOutput(command//SKC_" "//getPathVerbatimPosix(path), list, errmsg_def, exitstat = exitstat)
            if (failed) then
                if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg_def) ! LCOV_EXCL_LINE
                return ! LCOV_EXCL_LINE
            elseif (exitstat /= 0_IK) then
                index = reshape([integer(IK) ::], [0,0])
                list = SKC_""
                return
            end if
#define     UNESCAPE_REINDEX_LIST call setAsciiFromEscaped(list(index(1, i) : index(2, i)), endloc); index(2, i) = index(1, i) + endloc - 1_IK;
            call setSplit(index, list, NLC)
            lenList = size(index, 2, IK)
            ! Unescape the C-style escape sequences.
            if (showdir_def .and. showfile_def) then
                do i = 1_IK, lenList
                    UNESCAPE_REINDEX_LIST
                    !if (list(index(2, i) : index(2, i)) == shell%dirsep) index(2, i) = index(2, i) - 1_IK
                end do
            elseif (showfile_def) then
                do istart = 1_IK, lenList
                    if (list(index(2, i) : index(2, i)) /= shell%dirsep) exit
                end do
                do i = istart, lenList
                    UNESCAPE_REINDEX_LIST
                end do
                index = index(:, istart : lenList)
            elseif (showdir_def) then
                do istart = 1_IK, lenList
                    if (list(index(2, i) : index(2, i)) /= shell%dirsep) exit
                    UNESCAPE_REINDEX_LIST
                    !if (list(index(2, i) : index(2, i)) == shell%dirsep) index(2, i) = index(2, i) - 1_IK
                end do
                index = index(:, 1 : istart)
            else
                index = reshape([integer(IK) ::], [0,0])
                list = SKC_""
                return
            end if
#undef      UNESCAPE_REINDEX_LIST

        elseif (shell%is%windows .and. shell%is%cmd) then

            ! see: https://www.computerhope.com/dirhlp.htm
            !   /p                  Displays one screen of the listing at a time. To see the next screen, press any key.
            !   /q                  Displays file ownership information.
            !   /w                  Displays the listing in wide format, with as many as five file names or directory names on each line.
            !   /d                  Displays the listing in the same format as /w, but the files are sorted by column.
            !   /a[[:]<attributes>] Displays only the names of those directories and files with your specified attributes. If you don't use this parameter, the command displays the names of all files except hidden and system files. If you use this parameter without specifying any attributes, the command displays the names of all files, including hidden and system files. The list of possible attributes values are:
            !                       d - Directories
            !                       h - Hidden files
            !                       s - System files
            !                       l - Reparse points
            !                       r - Read-only files
            !                       a - Files ready for archiving
            !                       i - Not content indexed files
            !                       You can use any combination of these values, but do not separate your values using spaces.
            !                       Optionally you can use a colon (:) separator, or you can use a hyphen (-) as a prefix to mean, "not".
            !                       For example, using the -s attribute will not show the system files.
            !   /o[[:]<sortorder>]  Sorts the output according to sortorder, which can be any combination of the following values:
            !                       n - Alphabetically by name
            !                       e - Alphabetically by extension
            !                       g - Group directories first
            !                       s - By size, smallest first
            !                       d - By date/time, oldest first
            !                       Use the - prefix to reverse the sort order
            !                       Multiple values are processed in the order in which you list them. Do not separate multiple values with spaces, but you can optionally use a colon (:).
            !                       If sortorder is not specified, dir /o lists the directories alphabetically, followed by the files, which are also sorted alphabetically.
            !   /t[[:]<timefield>]  Specifies which time field to display or to use for sorting. The available timefield values are:
            !                       c - Creation
            !                       a - Last accessed
            !                       w - Last modified
            !   /s                  Lists every occurrence of the specified file name within the specified directory and all subdirectories.
            !   /b                  Displays a bare list of directories and files, with no additional information. The /b parameter overrides /w.
            !   /l                  Displays unsorted directory names and file names, using lowercase.
            !   /n                  Displays a long list format with file names on the far right of the screen.
            !   /x                  Displays the short names generated for non-8dot3 file names. The display is the same as the display for /n, but the short name is inserted before the long name.
            !   /c                  Displays the thousand separator in file sizes. This is the default behavior. Use /-c to hide separators.
            !   /4                  Displays years in four-digit format.
            !   /r                  Display alternate data streams of the file.
            !   /?                  Displays help at the command prompt.
            !
            !   \see
            !   https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/dir
            !
            command = SKC_"dir /b /o:g"
            if (showhidden_def) then
                command = command//SKC_" /a:h"
            else
                command = command//SKC_" /a:-h"
            end if
            if (present(sort)) then
                if (sort /= SK_"name") then
                    if (reversed_def) then
                        command = command//SKC_" /o:-n"
                    else
                        command = command//SKC_" /o:n"
                    end if
                elseif (sort == SK_"tmod") then
                    if (reversed_def) then
                        command = command//SKC_" /o:-d /t:w"
                    else
                        command = command//SKC_" /o:d /t:w"
                    end if
                elseif (sort == SK_"tacc") then
                    if (reversed_def) then
                        command = command//SKC_" /o:-d /t:a"
                    else
                        command = command//SKC_" /o:d /t:a"
                    end if
                elseif (sort == SK_"size") then
                    if (reversed_def) then
                        command = command//SKC_" /o:-s"
                    else
                        command = command//SKC_" /o:s"
                    end if
                elseif (sort == SK_"fext") then
                    if (reversed_def) then
                        command = command//SKC_" /o:-e"
                    else
                        command = command//SKC_" /o:e"
                    end if
                else
                    if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Unrecognized value for the `sort` input argument: "//sort
                    failed = .true._LK
                    return
                end if
            else
                if (reversed_def) then
                    command = command//SKC_" /o:-n"
                else
                    command = command//SKC_" /o:n"
                end if
            end if
            !if (len(exclude)>0 ) command = command // " | findstr /v /i " // exclude

            if (showdir_def) then
                if (isFailedGetOutput(command//SKC_" /a:d"//getPathVerbatimCMD(path), dirlist, errmsg_def, exitstat = exitstat)) then
                    dirIndex = reshape([integer(IK) ::], [0,0]) ! LCOV_EXCL_LINE
                    dirlist = SKC_"" ! LCOV_EXCL_LINE
                elseif (exitstat /= 0_IK) then
                    index = reshape([integer(IK) ::], [0,0])
                    list = SKC_""
                    return
                else
                    call setSplit(dirIndex, dirlist, NLC)
                    !do i = 1_IK, size(dirIndex,2,IK) - 1_IK
                    !    dirIndex(2, i) = dirIndex(2, i) + 1_IK
                    !    dirlist(dirIndex(2, i) : dirIndex(2, i)) = shell%dirsep
                    !end do
                    !dirIndex(2, i) = dirIndex(2, i) + 1_IK
                    !dirlist = dirlist//shell%dirsep
                end if
            end if

            if (showfile_def) then
                if (isFailedGetOutput(command//SKC_" /a:-d"//getPathVerbatimCMD(path), filelist, errmsg_def, exitstat = exitstat)) then
                    fileIndex = reshape([integer(IK) ::], [0,0]) ! LCOV_EXCL_LINE
                    filelist = SKC_"" ! LCOV_EXCL_LINE
                elseif (exitstat /= 0_IK) then
                    index = reshape([integer(IK) ::], [0,0])
                    list = SKC_""
                    return
                else
                    call setSplit(fileIndex, filelist, NLC)
                end if
            end if

            if (len(dirlist, IK) > 0_IK) then
                lenList = len(dirlist,IK) + 1_IK + len(filelist,IK) ! 1_IK takes care of the additional backslash for the last folder.
            else
                lenList = len(dirlist,IK) + len(filelist,IK)
            end if

            allocate(character(lenList,SKC) :: list)
            allocate(index(2,lenList))
            do i = 1_IK, size(dirIndex,2,IK)
                index(1, i) = dirIndex(1, i)
                index(2, i) = dirIndex(2, i) + 1_IK
                list(dirIndex(1, i) : dirIndex(2, i)) = dirlist(dirIndex(1, i) : dirIndex(2, i))
                list(index(2, i) : index(2, i)) = shell%dirsep
            end do

            istart = index(2, i)
            do i = 1_IK, size(fileIndex,2,IK)
                index(1,istart+i) = fileIndex(1, i)
                index(2,istart+i) = fileIndex(2, i)
                list(index(1,istart) + fileIndex(1, i) : index(1,istart) + fileIndex(2, i)) = filelist(fileIndex(1, i) : fileIndex(2, i))
            end do

        elseif (shell%is%windows .and. shell%is%powershell) then

            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": PowerShell is currently unsupported." ! LCOV_EXCL_LINE
            failed = .true._LK ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE

        else

            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Unrecognized runtime shell." ! LCOV_EXCL_LINE
            failed = .true._LK ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE

        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedCopy_ENABLED || isFailedMove_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*,SKC), parameter     :: DSP = DIR_SEP_POSIX
        character(*,SKC), parameter     :: DSW = DIR_SEP_WINDOWS
        character(*,SKC), parameter     :: DSPA = DIR_SEP_POSIX_ALL
        character(*,SKC), parameter     :: DSWA = DIR_SEP_WINDOWS_ALL
        character(*,SKC), parameter     :: LF = new_line(SKC_"a")
        character(:,SKC), allocatable   :: paths
        character(:,SKC), allocatable   :: dirname
        character(:,SKC), allocatable   :: command
        logical(LK)                     :: wait_def
        logical(LK)                     :: forced_def
        integer(IK)                     :: ntry_def
        integer(IK)                     :: itry
        type(shellis_type)              :: shellis
#if     isFailedMove_ENABLED
        character(*,SKC), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedMove()"
#elif   isFailedCopy_ENABLED
        character(*,SKC), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@isFailedCopy()"
        logical(LK)                     :: recursive_def
        if (present(recursive)) then
            recursive_def = recursive
        else
            recursive_def = .false._LK
        end if
#endif

        ERROR_STOP_IF(len(DSWA) /= 2, PROCEDURE_NAME//SK_": Internal library error occurred. The condition `len(DSWA) == 2` must hold.")

        ! The input source and destination cannnot be empty.

        failed = logical(from == SKC_"" .or. to == SKC_"", LK)
        if (failed) then
            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": The input path `from` is empty."
            return
        end if

        if (present(forced)) then
            forced_def = forced
        else
            forced_def = .false._LK
        end if

        shellis = shellis_type(failed)

        ! Define the copy/move command

        blockDefineCommand: if (shellis%posix .or. shellis%fish) then

            ! Create the destination nested directory if needed.

            if (isDir(from)) then
                ! `to` must be also a directory. If it does not exist, create it.
                if (isFile(to)) then
                    if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": The input destination `to` must be a directory when the source `from` is a directory. to: """//getStr(to)//SK_"""" ! LCOV_EXCL_LINE
                    failed = .true._LK ! LCOV_EXCL_LINE
                    return ! LCOV_EXCL_LINE
                elseif (.not. isDir(to)) then
                    dirname = getDirName(to, DSPA)
                    if (.not. isDir(dirname)) then ! it is a nonexistent dirname.
                        failed = isFailedMakeDir(dirname, errmsg = errmsg)
                        if (failed) then; if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg); return; end if ! LCOV_EXCL_LINE
                    end if
                end if
            else! `to` may still be a directory. If it is a directory, ensure it exists before copying the file.
                if (to(len(to):len(to)) == DSPA) then ! it must be a dir.
                    if (.not. isDir(to)) then ! it is a nonexistent (potentially nested) dir. Create it.
                        failed = isFailedMakeDir(to, errmsg = errmsg)
                        if (failed) then; if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg); return; end if ! LCOV_EXCL_LINE
                    end if
                else ! `to` must be a file. So make sure its (potentially nested) dirname exists.
                    dirname = getDirName(to, DSPA)
                    if (.not. isDir(dirname)) then ! it is a nonexistent dirname.
                        failed = isFailedMakeDir(dirname, errmsg = errmsg)
                        if (failed) then; if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg); return; end if ! LCOV_EXCL_LINE
                    end if
                end if
            end if

            ! quote the paths.

            if (shellis%powershell) then
                paths = getPathVerbatimPowerShell(from)//SKC_" "//getPathVerbatimPowerShell(to)
            else
                paths = getPathVerbatimPosix(from)//SKC_" "//getPathVerbatimPosix(to)
            end if

            ! copy/move path to destination.

#if         isFailedMove_ENABLED
            !   -n, --no-clobber    : do not overwrite an existing file.
            !   -f, --force         : if an existing destination file cannot be opened, remove it and try again (this option is ignored when the -n option is also used).
            if (forced_def) then
                command = SKC_"mv -f "//paths
            else
                command = SKC_"mv -n "//paths
            end if
#elif       isFailedCopy_ENABLED
            !   -a, --archive       : same as -dR --preserve=all.
            !   -d                  : same as --no-dereference --preserve=links
            !   -f, --force         : if an existing destination file cannot be opened, remove it and try again (this option is ignored when the -n option is also used).
            !   -n, --no-clobber    : do not overwrite an existing file.
            !   -R, -r, --recursive : copy directories recursively.
            if (forced_def .and. recursive_def) then
                command = SKC_"cp -arf "//paths
            elseif (recursive_def) then
                command = SKC_"cp -arn "//paths
            elseif (forced_def) then
                command = SKC_"cp -f "//paths
            else
                command = SKC_"cp "//paths
            end if
#endif

        elseif (shellis%windows) then blockDefineCommand ! either CMD or PowerShell: robocopy automatically generates new directories if needed.

            ! LCOV_EXCL_START
            if (isDir(from)) then ! use robocopy.

                ! To get the linux `cp` behavior on Windows when `to` is an existing directory, we should create a new directory in `to` with the basename of `from`.
                if (isDir(to)) then
                    paths = getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to//DSW//getBaseName(from, DSWA))
                else
                    paths = getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to)
                end if

                !   Robocopy flags:
                !   /e          :   Copies subdirectories. This option automatically includes empty directories.
                !   /r:<n>      :   Specifies the number of retries on failed copies. The default value of n is 1,000,000 (one million retries).
                !   /w:<n>      :   Specifies the wait time between retries, in seconds. The default value of n is 30 (wait time 30 seconds).
                !   /xc         :   Excludes changed files.
                !   /xn         :   Excludes newer files.
                !   /xo         :   Excludes older files.
                !   /xc /xn /xo :   Does not overwrite existing files.
                !   /is         :   Includes the same files. Same files are identical in name, size, times, and all attributes.
                !   /move       :   Moves files and directories, and deletes them from the source after they are copied.
#if             isFailedMove_ENABLED
                if (forced_def) then
                    command = SKC_"robocopy "//paths//SKC_" /r:1 /w:1 /move /is         /e"//SKC_" > nul"
                else
                    command = SKC_"robocopy "//paths//SKC_" /r:1 /w:1 /move /xc /xn /xo /e"//SKC_" > nul"
                end if
#elif           isFailedCopy_ENABLED
                if (forced_def .and. recursive_def) then
                    command = SKC_"robocopy "//paths//SKC_" /r:1 /w:1 /is         /e"//SKC_" > nul"
                elseif (recursive_def) then
                    command = SKC_"robocopy "//paths//SKC_" /r:1 /w:1 /xc /xn /xo /e"//SKC_" > nul"
                elseif (forced_def) then
                    command = SKC_"robocopy "//paths//SKC_" /r:1 /w:1 /is           "//SKC_" > nul"
                else
                    command = SKC_"robocopy "//paths//SKC_" /r:1 /w:1 /xc /xn /xo   "//SKC_" > nul"
                end if
#endif

            else! `from` must be a file. Use xcopy/move, which can handle files with nested destinations.

                if (.not. forced_def) then
                    if (isFile(to)) then
                        return ! Attempting to overwrite a single existing file. Nothing to do. Return.
                    elseif (isDir(to)) then
                        if (isFile(to//DSW//getBaseName(from, DSWA))) return ! Attempting to overwrite a single existing file whose directory is `to`. Nothing to do. Return.
                    end if
                end if

                ! Rest assured; `to` is not an existing file. However, `to` may still be intended as a (potentially non-existing) directory if it ends with a directory separator.
                ! If so, pipe D to xcopy/move indicate it is a directory.

                ! copy or xcopy flags:
                ! /e Copies all (even empty) subdirectories.
                ! /q Suppresses the display of xcopy messages.
                ! /y Suppresses prompting to confirm that you want to overwrite an existing destination file (as if `forced = .true.`).
                ! /s Copies directories and subdirectories, unless they are empty. If one omits /s, xcopy works within a single directory.
                ! command = SKC_"xcopy /e /q /s /Y "//getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to)//SKC_" > nul"

                itry = len(to, IK)
                if (index(DSWA, to(itry:itry), kind = IK) > 0_IK) then ! it must be a dir.
                    ! \warning
                    ! `xcopy` cannot automatically create (nested) folders that end with a forward slash (no problem with backward slash.
                    ! Therefore, a safer and faster way of creating the directory is to let xcopy do it, while ignoring the last character (which is a directory separator).
                    ! `echo D` pipes the response `directory` to the question of how to interpret `to` if there is ambiquity (e.g., if `to` does not end with directory separator).
#if                 isFailedMove_ENABLED
                    ! move can neither realize `to` is a directory (on powershell) nor create a (potentially nested) `to` folder (in either cmd or powershell). So create it.
                    failed = isFailedMakeDir(to, errmsg = errmsg)
                    if (failed) then; if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg); return; end if ! LCOV_EXCL_LINE
#elif               isFailedCopy_ENABLED
                    command = SKC_"echo D | xcopy /q /y "//getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to(1:itry-1))//SKC_" > nul"
                else
                    ! `to` must be a file. So no problemino.
                    ! `echo F` pipes the response `file` to the question of how to interpret `to` if there is ambiquity (e.g., if `to` does not end with directory separator).
                    command = SKC_"echo F | xcopy /q /y "//getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to)//SKC_" > nul"
#endif
                end if

#if             isFailedMove_ENABLED
                if (shellis%cmd) then
                    command = SKC_"move /y "//getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to)//SKC_" > nul"
                elseif (shellis%powershell) then
                    command = SKC_"move "//getPathVerbatimCMD(from)//SKC_" "//getPathVerbatimCMD(to)//SKC_" -Force > nul"
                end if
#endif
            end if
            ! LCOV_EXCL_STOP

        else blockDefineCommand

            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Failed to fetch the runtime shell type."//LF//trim(errmsg) ! LCOV_EXCL_LINE
            failed = .true._LK ! LCOV_EXCL_LINE
            return ! LCOV_EXCL_LINE

        end if blockDefineCommand

        ! Attempt repeatedly to copy/move the file.

        if (present(ntry)) then
            CHECK_ASSERTION(__LINE__, 0_IK < ntry, PROCEDURE_NAME//SK_": The condition `0 < ntry` must hold. ntry = "//getStr(ntry))
            ntry_def = ntry
        else
            ntry_def = 1_IK
        end if

        if (present(wait)) then
            wait_def = wait
        else
            wait_def = .true._LK
        end if

        do itry = 1_IK, ntry_def
            failed = isFailedExec(command, wait = wait_def, cmdmsg = errmsg)
            !failed = wait_def .and. .not. isExtant(to)
            if (.not. failed) exit
        end do

        if (failed .and. present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Failed to accomplish task file after "//getStr(ntry_def)//SK_" attempts."//LF//trim(errmsg)

        !   \warning do NOT remove.
        !   \bug gfortran with heap memory allocations fails to automatically deallocate on exit.
        deallocate(command)

        !%%%%%%%%%%%%%%%%%%%%%
#elif   isFailedRemove_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//"@isFailedRemove()"
        character(*,SKC), parameter :: LF = new_line(SKC_"a")
        logical(LK) :: forced_def, recursive_def, wait_def
        character(:,SKC), allocatable :: command, cmdfull
        character(511, SK) :: errmsg_def
        integer(IK) :: i, itry, ntry_def
        integer(IK) :: iostat

        errmsg_def = ""
        failed = .false.
        ntry_def = 1_IK
        wait_def = .true._LK
        forced_def = .false._LK
        recursive_def = .false._LK
        if (present(ntry)) ntry_def = ntry
        if (present(wait)) wait_def = wait
        if (present(forced)) forced_def = forced
        if (present(recursive)) recursive_def = recursive
        CHECK_ASSERTION(__LINE__, 0_IK < ntry_def, PROCEDURE_NAME//SK_": The condition `0 < ntry` must hold. ntry = "//getStr(ntry_def))

        ! non-existence is okay only if `forced` is `.true.`.

        if (isFile(path)) then
            fileRemoval: block
                logical(LK) :: opened
                integer(IK) :: unit
                inquire(file = path, opened = opened, number = unit, iostat = iostat, iomsg = errmsg_def)
                if (iostat /= 0_IK) exit fileRemoval
                if (.not. opened) then
                    open(newunit = unit, file = path, status = "replace", iostat = iostat, iomsg = errmsg_def SHARED)
                    if (iostat /= 0_IK) exit fileRemoval
                end if
                do itry = 1, ntry_def
                    close(unit, status = "delete", iostat = iostat, iomsg = errmsg_def)
                    if (iostat /= 0_IK) then
                        errmsg_def = SK_"The `close(unit,status='delete')` statement failed. "//trim(errmsg_def)
                        exit fileRemoval
                    end if
                    inquire(file = path, exist = failed, iostat = iostat, iomsg = errmsg_def)
                    if (iostat /= 0_IK) exit fileRemoval
                    failed = wait_def .and. failed
                    if (.not. failed) return
                end do
                errmsg_def = SK_": Failed to accomplish task file after "//getStr(ntry_def)//SK_" attempts."
            end block fileRemoval
            failed = .true._LK
            if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg_def)
            return
        end if

        ! Remove possibly a directory or a pattern.
        ! Note that some commands like Windows CMD `rd` do not digest wildcards.
        ! As such a list of paths corresponding to the specified wildcards must be created first.

        recursiveRemoval: block

            type(shellis_type) :: shellis
            type(css_type), allocatable :: list(:), verbatim(:)

            if (isDir(path)) then
                list = [css_type(path)]
            else
                failed = isFailedGlob(path, list, errmsg_def)
                if (failed) exit recursiveRemoval
                if (size(list, 1, IK) == 0_IK) then
                    failed = .not. forced_def
                    if (failed) then
                        errmsg_def = SK_"The specified path does not match any file or directory. Set `forced = .true.` to gracefully ignore non-existing paths."
                        exit recursiveRemoval
                    else
                        return
                    end if
                end if
            end if

            ! Removing non-file paths requires `recursive` option to be `.true.`.

            if (0_IK < size(list, 1, IK) .and. .not. recursive_def) then
                if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": Removing a directory or a set of paths matching a pattern requires setting the input argument `recursive = .true.`."
                failed = .true._LK
                return
            end if

            shellis = shellis_type(failed, errmsg_def)
            call setResized(verbatim, size(list, 1, IK))

            blockDefineCommand: if (shellis%posix .or. shellis%fish) then

                do i = 1, size(list, 1, IK)
                    if (shellis%powershell) then
                        verbatim(i)%val = getPathVerbatimPowerShell(list(i)%val)
                    else
                        verbatim(i)%val = getPathVerbatimPosix(list(i)%val)
                    end if
                end do
                !   -f, --force         :   ignore nonexistent files and arguments, never prompt.
                !   -r, -R, --recursive :   remove directories and their contents recursively.
                command = SKC_"rm "
                if (forced_def) command = command//SKC_"-f "
                if (recursive_def) command = command//SKC_"-r "

            elseif (shellis%windows .and. shellis%cmd) then blockDefineCommand

                do i = 1, size(list, 1, IK)
                    verbatim(i)%val = getPathVerbatimCMD(list(i)%val)
                end do
                !   /Q  :   Quiet mode, do not ask if ok to remove a directory tree with /S.
                !   /S  :   Removes all directories and files in the specified directory in addition to the directory itself. Used to remove a directory tree.
                command = SKC_"rd "
                if (forced_def) command = command//SKC_"/Q "
                if (recursive_def) command = command//SKC_"/S "

            elseif (shellis%windows .and. shellis%powershell) then blockDefineCommand

                do i = 1, size(list, 1, IK)
                    verbatim(i)%val = getPathVerbatimPowerShell(list(i)%val)
                end do
                !   -Force      :   Forces the cmdlet to remove items that can't otherwise be changed, such as hidden or read-only files or read-only aliases or variables.
                !   -Recurse    :   Indicates that this cmdlet deletes the items in the specified locations and in all child items of the locations.
                !                   The Recurse parameter might not delete all subfolders or all child items. This is a known issue.
                command = SKC_"Remove-Item "
                if (forced_def) command = command//SKC_"-Force "
                if (recursive_def) command = command//SKC_"-Recurse "

            else blockDefineCommand

                errmsg_def = SK_": Failed to fetch the runtime shell type."//LF//trim(errmsg) ! LCOV_EXCL_LINE
                exit recursiveRemoval

            end if blockDefineCommand

            ! Remove paths recursively.

            loopOverPath: do i = 1, size(list, 1, IK)
                cmdfull = command//verbatim(i)%val
                do itry = 1, ntry_def
                    failed = isFailedExec(cmdfull, wait = wait_def, cmdmsg = errmsg_def)
                    inquire(file = list(i)%val, exist = failed, iostat = iostat, iomsg = errmsg_def)
                    if (iostat /= 0_IK) exit recursiveRemoval
                    failed = wait_def .and. failed
                    if (.not. failed) cycle loopOverPath
                end do
                errmsg_def = SK_": Failed to accomplish remove path '"//list(i)%val//SK_"' after "//getStr(ntry_def)//SK_" attempts."
                exit recursiveRemoval
            end do loopOverPath

            return

        end block recursiveRemoval

        if (present(errmsg)) errmsg = PROCEDURE_NAME//SK_": "//trim(errmsg_def)
        failed = .true._LK

        !%%%%%%%%%%%%%%%%%%%
#elif   getPathMatch_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        use pm_err, only: getFine
        integer(IK) :: lenList
        logical(LK) :: failed_def
        character(:,SKC), allocatable :: key_def, inc_def, sep_def, paths
        if (present(failed)) failed = .false._LK
        if (present(key)) then
            key_def = key
        else
            key_def = ""
        endif
        if (present(inc)) then
            inc_def = inc
        else
            inc_def = ""
        endif
        if (present(sep)) then
            sep_def = sep
        else
            sep_def = getPathSep(failed = failed_def, errmsg = errmsg)
            if (failed_def) then
                if (present(failed)) then
                    if (present(errmsg)) errmsg = getFine(__FILE__, __LINE__)//MODULE_NAME//"@getPathMatch(): Failed to infer the path separator character. "//trim(errmsg)
                    call setResized(list, 0_IK)
                    failed = failed_def
                    return
                else
                    error stop getFine(__FILE__, __LINE__)//MODULE_NAME//"@getPathMatch(): Failed to infer the path separator character. "//trim(errmsg)
                end if
            end if
        endif
        lenList = 4095_IK
        do
            paths = repeat(" ", lenList)
            call setPathMatch(key_def, inc_def, sep_def, paths, lenList)
            if (len(paths, IK) < abs(lenList)) then
                lenList = abs(lenList)
                cycle
            elseif (0_IK <= lenList) then
                call setSplit(list, paths(1 : lenList), sep_def)
                exit
            elseif (present(failed)) then
                if (present(errmsg)) errmsg = getFine(__FILE__, __LINE__)//MODULE_NAME//"@getPathMatch(): Failed to detect paths. "//trim(errmsg)
                call setResized(list, 0_IK)
                failed = .true._LK
                return
            else
                error stop getFine(__FILE__, __LINE__)//MODULE_NAME//"@getPathMatch(): Failed to infer the path separator character. "//trim(errmsg)
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%
#elif   setPathMatch_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ipath, lenkey, leninc, item, lenout, lenListInput, istart
        logical(LK) :: isall, found, failed
        character(255,SKC) :: errmsg
        character(1,SKC) :: pathsep
        !integer(IK), allocatable :: cindex(:,:)
        !integer(IK), allocatable :: kindex(:,:)
        type(css_type), allocatable :: csskey(:)
        type(css_type), allocatable :: inclist(:)
        type(css_type), allocatable :: csspath(:)
        character(:,SKC), allocatable :: contents, lower, inclow
        CHECK_ASSERTION(__LINE__, 0_IK < lenList, SK_"@setPathMatch(): The condition `0 < lenList` must hold. lenList = "//getStr(lenList))
        CHECK_ASSERTION(__LINE__, 0_IK < len_trim(sep), SK_"@setPathMatch(): The condition `0 < len_trim(sep)` must hold. len_trim(sep) = "//getStr(len_trim(sep)))
        lenListInput = lenList ! Keep a copy of the input value.
        errmsg = ""

        ! First search the PATH env variable.

        success: block
            failed = isFailedGetEnvVar("PATH", contents, errmsg)
            if (failed) exit success
            pathsep = getPathSep(failed, errmsg)
            if (.not. failed) then
#if             FORTRAN_ENABLED
                leninc = len(inc, IK)
                lenkey = len(key, IK)
#else
                lenkey = 1_IK
                do
                    if (key(lenkey : lenkey) == c_null_char) exit
                    lenkey = lenkey + 1_IK
                end do
                lenkey = lenkey - 1_IK
                leninc = 1_IK
                do
                    if (inc(leninc : leninc) == c_null_char) exit
                    leninc = leninc + 1_IK
                end do
                leninc = leninc - 1_IK
#endif
                lenList = 0_IK
                if (0_IK < leninc) inclow = getStrLower(inc(1 : leninc))
                call setSplit(csskey, getStrLower(key(1 : lenkey)), sep)
                isall = 0_IK == size(csskey, 1, IK)
                if (.not. isall) then
                    call setSplit(csspath, contents, pathsep)
                    isall = 0_IK == size(csspath, 1, IK)
                    do ipath = 1, size(csspath, 1, IK)
                        lower = getStrLower(csspath(ipath)%val)
                        found = .true._LK
                        do item = 1, size(csskey, 1, IK)
                            found = 0 < index(lower, csskey(item)%val)
                            if (.not. found) exit
                        end do
                        if (found) then
                            ! Add the path only if the `inc` is also present.
                            found = 0_IK == leninc
                            if (.not. found) then
                                failed = isFailedList(csspath(ipath)%val, inclist, errmsg = errmsg)
                                if (failed) exit success
                                do item = 1, size(inclist, 1, IK)
                                    found = 0 < index(getStrLower(inclist(item)%val), inclow)
                                    if (found) exit
                                end do
                            end if
                            if (found) then
                                ! Add the separator.
                                if (0_IK < lenList) then
                                    if (lenList + 1 <= lenListInput) list(lenList + 1 : lenList + 1) = sep
                                    lenList = lenList + 1_IK
                                end if
                                ! Add the path.
                                if (lenList + len(lower, IK) <= lenListInput) list(lenList + 1 : lenList + len(lower, IK)) = csspath(ipath)%val
                                lenList = lenList + len(lower, IK)
                            end if
                        end if
                    end do
                end if
                if (isall) then
                    call setReplaced(contents, pathsep, sep)
                    lenout = len_trim(contents)
                    lenList = min(lenout, lenListInput)
                    list(1 : lenList) = contents(1 : lenList)
                    lenList = merge(lenout, lenList, lenListInput < lenout)
                end if
            end if

            if (isDarwin()) then
                failed = isFailedGetOutput(SK_"system_profiler SPApplicationsDataType", contents, errmsg)!, exitstat = exitstat)
                if (failed) exit success
                ! The following is the typical output of the command:
                !
                !    MATLAB_R2019a:
                !
                !          Version: R2019a (9.6.0)
                !          Obtained from: Unknown
                !          Last Modified: 6/18/21, 2:15 AM
                !          Location: /Applications/MATLAB_R2019a.inc
                !          Kind: 64-bit
                !
                ! We search for all instances of lines matching "Location: " and the input keys.
                call setSplit(csspath, contents, new_line(SKC_""))
                do ipath = 1, size(csspath, 1, IK)
                    istart = index(csspath(ipath)%val, SK_"Location: ", kind = IK)
                    if (0 < istart) then
                        istart = istart + 10_IK
                        lower = getStrLower(csspath(ipath)%val(istart :))
                        found = .true._LK
                        do item = 1, size(csskey, 1, IK)
                            found = 0 < index(lower, csskey(item)%val)
                            if (.not. found) exit
                        end do
                        if (found) then
                            ! Add the path only if the `inc` is also present.
                            found = 0_IK == leninc
                            if (.not. found) then
                                failed = isFailedList(csspath(ipath)%val(istart :), inclist, errmsg = errmsg)
                                if (failed) exit success
                                do item = 1, size(inclist, 1, IK)
                                    found = 0 < index(getStrLower(inclist(item)%val), inclow)
                                    if (found) exit
                                end do
                            end if
                            if (found) then
                                ! Add the separator.
                                if (0_IK < lenList) then
                                    if (lenList + 1 <= lenListInput) list(lenList + 1 : lenList + 1) = sep
                                    lenList = lenList + 1_IK
                                end if
                                ! Add the path.
                                if (lenList + len(lower, IK) <= lenListInput) list(lenList + 1 : lenList + len(lower, IK)) = csspath(ipath)%val(istart :)
                                lenList = lenList + len(lower, IK)
                            end if
                        end if
                    end if
                end do
            end if

        end block success

        !if (isWindows() .or. isLinux()) then
        if (failed) then
            lenout = len_trim(errmsg)
            lenList = min(lenout, lenListInput)
            list(1 : lenList) = errmsg(1 : lenList)
            list(lenList + 1 : lenListInput) = ""
            lenList = -merge(lenout, lenList, lenListInput < lenout)
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  SHARED