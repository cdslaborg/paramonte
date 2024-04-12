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
!>  This file contains procedure implementations of [pm_fftpack](@ref pm_fftpack).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the compiler specs.
#if     INTEL_ENABLED && WINDOWS_ENABLED
#define INTEL_SHARED_FILE, SHARED
#else
#define INTEL_SHARED_FILE
#endif
        ! Define the resource.
#if     File_ENABLED
#define ITEM file
#elif   Unit_ENABLED || Unit_ENABLED
#define ITEM unit
#elif   isOpen_ENABLED || getAction_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the runtime error check.
#if     (setContentsTo_ENABLED || setContentsFrom_ENABLED) && CII_ENABLED
#define RETURN_IF_FAILED if (iostat /= 0_IK) return ! LCOV_EXCL_LINE
#define IOSTAT_IOMSG , iostat = iostat, iomsg = iomsg
#elif   (setContentsTo_ENABLED || setContentsFrom_ENABLED) && CDD_ENABLED
#define RETURN_IF_FAILED
#define IOSTAT_IOMSG
#elif   (getErrTableRead_ENABLED || getErrTableWrite_ENABLED) && File_ENABLED
#define RETURN_IF_FAILED(LINE) if (err /= 0_IK) then; if (present(iomsg)) iomsg = getStr(LINE)//SK_": "//iomsg_def; close(unit); return; end if ! LCOV_EXCL_LINE
#elif   (getErrTableRead_ENABLED || getErrTableWrite_ENABLED) && Unit_ENABLED
#define RETURN_IF_FAILED(LINE) if (err /= 0_IK) then; if (present(iomsg)) iomsg = getStr(LINE)//SK_": "//iomsg_def; return; end if ! LCOV_EXCL_LINE
#elif   (setContentsTo_ENABLED || setContentsFrom_ENABLED || getErrTableRead_ENABLED || getErrTableWrite_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Define the error check.
#define SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) \
if (present(iostat)) then; iostat = iostat_def; if (present(iomsg)) iomsg = iomsg_def; elseif (iostat_def /= 0_IK) then; error stop SK_"FATAL RUNTIME ERROR: "//trim(adjustl(iomsg_def)); end if;

        !%%%%%%%%%%%%%
#if     isOpen_ENABLED
        !%%%%%%%%%%%%%

        inquire(ITEM = ITEM, opened = opened)

        !%%%%%%%%%%%%%%%%
#elif   getAction_ENABLED
        !%%%%%%%%%%%%%%%%

        inquire(ITEM = ITEM, action = action)

        !%%%%%%%%%%%%%%%%%%%%%
#elif   constructField_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        if (present(string)) field%string = string
        if (present(integer)) field%integer = integer
        if (present(logical)) field%logical = logical
        if (present(complex)) field%complex = complex
        if (present(real)) field%real = real

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getCountRecord_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

#if     File_ENABLED
        character(LEN_IOMSG)    :: iomsg_def
        integer(IK)             :: unit, iostat_def
        CHECK_ASSERTION(__LINE__, .not. isOpen(file), SK_"@getCountRecord(): The condition `.not. isOpen(file)` must hold. file = "//getStr(file))
        open( iomsg = iomsg_def & ! LCOV_EXCL_LINE
            , iostat = iostat_def & ! LCOV_EXCL_LINE
            , newunit = unit & ! LCOV_EXCL_LINE
            , position = "rewind" & ! LCOV_EXCL_LINE
            , access = "sequential" & ! LCOV_EXCL_LINE
            , action = "read" & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
            , file = file & ! LCOV_EXCL_LINE
            INTEL_SHARED_FILE)
        SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) ! LCOV_EXCL_LINE
        if (iostat_def /= 0_IK) return ! LCOV_EXCL_LINE
#elif   Unit_ENABLED
        rewind(unit)
        CHECK_ASSERTION(__LINE__, isOpen(unit), SK_"@getCountRecord(): The condition `isOpen(unit)` must hold. unit = "//getStr(unit))
#else
#error  "Unrecognized interface."
#endif
        nrecord = getCountRecordLeft(unit, isCountable, iostat = iostat, iomsg = iomsg) ! Count the file records.
#if     File_ENABLED
        ! Close / delete the file.
        call setFileClosed(unit, del, iostat, iomsg)
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCountRecordLeft_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%

        character(:, SK), allocatable   :: record
        character(LEN_IOMSG)            :: iomsg_def
        integer(IK)                     :: iostat_def
        integer(IK)                     :: ub

        ! Count the file records.

        nrecord = 0_IK
        if (present(isCountable)) then
            do
                call setRecordFrom(unit = unit, record = record, iostat = iostat_def, iomsg = iomsg_def, ub = ub)
                if (iostat_def == 0_IK) then
                    if (isCountable(record(1:ub))) nrecord = nrecord + 1_IK
                elseif (iostat_def == iostat_end) then
                    deallocate(record)
                    exit
                else
                    SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) ! LCOV_EXCL_LINE
                    return ! LCOV_EXCL_LINE
                end if
            end do
        else
            do
                read(unit, "(A)", iostat = iostat_def, iomsg = iomsg_def)
                if (iostat_def == 0_IK) then
                    nrecord = nrecord + 1_IK
                elseif (iostat_def == iostat_end) then
                    exit
                else
                    SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) ! LCOV_EXCL_LINE
                    return ! LCOV_EXCL_LINE
                end if
            end do
        end if

        if (present(reset)) then
            if (reset) then
                do ub = 1_IK, nrecord + 1_IK
                    backspace(unit, iostat = iostat_def)
                    if (iostat_def == 0_IK) cycle
                    iomsg_def = MODULE_NAME//SK_"@getCountRecordLeft(): Failed to backspace record."
                    SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) ! LCOV_EXCL_LINE
                    return ! LCOV_EXCL_LINE
                end do
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getContentsFrom_ENABLED && Unit_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setContentsFrom(unit = unit, contents = contents, del = del)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getContentsFrom_ENABLED && File_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setContentsFrom(file = file, contents = contents, del = del)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setContentsFrom_ENABLED && Unit_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*,SKC)    , parameter     :: LF = new_line(SKC_"a")
        integer(IK)         , parameter     :: LENLF = len(LF, IK)
        integer(IK)                         :: lb, ub, pos, size, recl
        character(10, SK)                   :: access
#if     CDD_ENABLED
#define CATCH_ERR_IF_FAILED(FINE) \
if (iostat /= 0_IK) error stop FINE//MODULE_NAME//SK_"@setContentsFrom(): "//trim(iomsg) ! LCOV_EXCL_LINE
        integer(IK)                         :: iostat
        character(LEN_IOMSG, SK)            :: iomsg
#elif   CII_ENABLED
#define CATCH_ERR_IF_FAILED(FINE) \
if (iostat /= 0_IK) return ! LCOV_EXCL_LINE
#else
#error  "Unrecognized interface."
#endif

        CHECK_ASSERTION(__LINE__, isOpen(unit), SK_"@setContentsFrom(): The condition `isOpen(unit)` must hold. unit = "//getStr(unit))
        CHECK_ASSERTION(__LINE__, index(getAction(unit), "READ") > 0, SK_"@setContentsFrom(): The condition `index(getAction(unit), ""READ"") > 0` must hold. unit = "//getStr(unit))

        inquire(unit = unit, access = access IOSTAT_IOMSG)
        RETURN_IF_FAILED

        if (access == SK_"STREAM") then
            inquire(unit = unit, pos = pos, size = size IOSTAT_IOMSG)
            RETURN_IF_FAILED
            allocate(character(size,SKC) :: contents)
            read(unit, pos = pos IOSTAT_IOMSG) contents
        elseif (access == SK_"DIRECT") then
            inquire(unit = unit, nextrec = pos, recl = recl IOSTAT_IOMSG)
            RETURN_IF_FAILED
            recl = recl + LENLF
            allocate(character(recl,SKC) :: contents)
            size = len(contents, IK)
            lb = 1_IK
            ub = recl
            do
                read(unit, rec = pos IOSTAT_IOMSG) contents(lb : ub - LENLF)
                if (iostat == iostat_end) exit
                CATCH_ERR_IF_FAILED(getFine(__FILE__, __LINE__))
                contents(lb + recl : ub) = LF
                lb = ub + 1_IK
                ub = ub + recl
                pos = pos + 1_IK
                if (size < ub) call setResized(contents, ub)
            end do
            contents = contents(1 : ub - LENLF)
#if         CII_ENABLED
            iostat = 0_IK
#endif
        else!if (access == SK_"SEQUENTIAL") then or it could be "UNDEFINED" if not set explicitly in gfortran 13.
            !error stop MODULE_NAME//SK_"@setContentsFromUnit(): An impossible Internal library error detected. The access attribute of the input `unit` is unrecognized. access="//access ! LCOV_EXCL_LINE
            allocate(character(LEN_IOMSG,SKC) :: contents)
            lb = 1_IK
            do
                call setRecordFrom(unit, contents, iostat, iomsg, lb = lb, ub = ub, linefed = .true._LK)
                if (iostat == iostat_end) exit
                CATCH_ERR_IF_FAILED(getFine(__FILE__, __LINE__))
                lb = ub + 1_IK
            end do
            contents = contents(1 : ub - LENLF)
#if         CII_ENABLED
            iostat = 0_IK
#endif
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setContentsFrom_ENABLED && File_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: unit
        integer(IK) :: lenContents
        CHECK_ASSERTION(__LINE__, .not. isOpen(file), SK_"@setContentsFrom(): The condition `.not. isOpen(file)` must hold. file = "//getStr(file))
        open( file = file & ! LCOV_EXCL_LINE
            , newunit = unit & ! LCOV_EXCL_LINE
            , form = "unformatted" & ! LCOV_EXCL_LINE
            , position = "rewind" & ! LCOV_EXCL_LINE
            , access = "stream" & ! LCOV_EXCL_LINE
            , action = "read" & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
            IOSTAT_IOMSG INTEL_SHARED_FILE)
        RETURN_IF_FAILED

        ! Inquire the file size in bytes.

        inquire(unit = unit, size = lenContents IOSTAT_IOMSG)
        RETURN_IF_FAILED

        ! Read the file contents as a string.

        allocate(character(lenContents, SK) :: contents)
        read(unit IOSTAT_IOMSG) contents
        RETURN_IF_FAILED

        ! Close/delete the file.

        if (present(del)) then
            if (del) then
                close(unit, status = "delete" IOSTAT_IOMSG)
                return
            end if
        end if
        close(unit IOSTAT_IOMSG)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setContentsTo_ENABLED && File_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: unit
        CHECK_ASSERTION(__LINE__, .not. isOpen(file), SK_"@setContentsFrom(): The condition `.not. isOpen(file)` must hold. file = "//getStr(file))
        open( file = file & ! LCOV_EXCL_LINE
            , newunit = unit & ! LCOV_EXCL_LINE
            , position = "rewind" & ! LCOV_EXCL_LINE
            , form = "unformatted" & ! LCOV_EXCL_LINE
            , access = "stream" & ! LCOV_EXCL_LINE
            , action = "write" & ! LCOV_EXCL_LINE
            , status = "replace" & ! LCOV_EXCL_LINE
            IOSTAT_IOMSG INTEL_SHARED_FILE)
        RETURN_IF_FAILED

        ! Write the contents to the file as a string.

        !write(unit, "(A)" IOSTAT_IOMSG) contents
        write(unit IOSTAT_IOMSG) contents
        RETURN_IF_FAILED
        close(unit IOSTAT_IOMSG)

        !%%%%%%%%%%%%%%%%%%%%
#elif   getRecordFrom_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        character(LEN_IOMSG, SK) :: iomsg_def
        integer(IK) :: iostat_def
        call setRecordFrom(unit, record, iostat_def, iomsg_def, linefed = linefed)
        if (present(iostat)) iostat = iostat_def
        if (present(iomsg)) iomsg = iomsg_def

        !%%%%%%%%%%%%%%%%%%%%
#elif   setRecordFrom_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: size, lb_def, lenRecord
#if     UR_ENABLED
        integer(IK) :: iostat
        character(LEN_IOMSG, SK) :: iomsg
#elif   !URII_ENABLED
#error  "Unrecognized interface."
#endif
        character(*,SKC), parameter :: LF = new_line(SKC_"a")
        integer(IK), parameter :: LENLF = len(LF, IK)
        integer(IK) :: nlen

        nlen = 0_IK
        if (present(linefed)) then
            if (linefed) nlen = LENLF
        end if

        if (present(lb)) then
            CHECK_ASSERTION(__LINE__, 0_IK < lb, SK_"@setRecordFrom(): The condition `0_IK < lb` must hold. lb = "//getStr(lb))
            lb_def = lb - 1_IK
        else
            lb_def = 0_IK
        end if

        if (allocated(record)) then
            lenRecord = len(record, IK) - nlen
            if (lenRecord <= lb_def) then
                deallocate(record)
                lenRecord = lb_def + LEN_RECORD
                allocate(character(lenRecord + nlen, SK) :: record)
            end if
        else
            lenRecord = lb_def + LEN_RECORD
            allocate(character(lenRecord + nlen, SK) :: record)
        end if

        do
            read(unit, "(a)", advance = "no", size = size, iostat = iostat, iomsg = iomsg) record(lb_def + 1_IK : lenRecord)
            if (iostat == iostat_eor) then ! Record reading is complete.
#if             URII_ENABLED
                iostat = 0_IK
                iomsg = getLine(__LINE__)//MODULE_NAME//SK_"@setRecordFrom(): "//trim(iomsg)
#endif
                lb_def = lb_def + size + nlen
                if (present(ub)) then
                    ub = lb_def
                else
                    call setResized(record, lb_def)
                end if
                if (nlen > 0_IK) record(lb_def - nlen + 1_IK : lb_def) = LF
                return
            elseif (iostat == 0_IK) then ! There is still record to read.
                lb_def = lb_def + size
                CHECK_ASSERTION(__LINE__, lb_def == lenRecord, SK_"@setRecordFrom(): Internal library error detected. The condition `lb_def == lenRecord` must hold. lb_def, lenRecord = "//getStr([lb_def, lenRecord]))
                lenRecord = lenRecord + lenRecord
                call setResized(record, lenRecord + nlen)
                cycle
            else
#if             UR_ENABLED
                error stop getLine(__LINE__)//MODULE_NAME//SK_"@setRecordFrom(): "//trim(iomsg)
#endif
                return
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%%
#elif   setFileClosed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        character(LEN_IOMSG)    :: iomsg_def
        integer(IK)             :: iostat_def, i
        if (present(del)) then
            if (del) then ! Attempt to delete the file repeatedly. This is important on windows systems as the file often remains locked.
                do i = 1_IK, 100_IK
                    close(unit, status = "delete", iostat = iostat_def, iomsg = iomsg_def)
                    if (iostat_def == 0_IK) then
                        if (present(iostat)) iostat = iostat_def
                        return
                    end if
                end do
                ! All attempts at closing and deleting the file failed.
                SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) ! LCOV_EXCL_LINE
                if (iostat_def /= 0_IK) return ! LCOV_EXCL_LINE
            end if
        end if
        close(unit, iostat = iostat_def, iomsg = iomsg_def)
        SET_STAT_IO(iostat_def, iomsg_def, iostat, iomsg) ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getErrTableRead_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     File_ENABLED
        integer(IK) :: unit
#endif
        character(LEN_IOMSG, SK) :: iomsg_def
        character(*, SK), parameter :: HT = achar(9, SK)
        character(*, SK), parameter :: LF = achar(10, SK)
        character(*, SK), parameter :: CR = achar(13, SK)
        integer(IK), parameter :: RINIT = 127_IK
        integer(IK) :: irow, nrow
#if     !(D1_ENABLED && NO_ENABLED)
        integer(IK), allocatable :: sepLoc(:)
        character(:, SK), allocatable :: record
        integer(IK) :: icol, ncol, lenrec, lenLoc, sepLen, nsep
#endif
#if     CK_ENABLED && !(D1_ENABLED && NO_ENABLED)
        real(CKC), allocatable :: field(:)
#endif
        ! Define the transposition rules.
#if     D2_ENABLED && NO_ENABLED
#define GET_INDEX(I,J) I,J
#elif   D2_ENABLED && TO_ENABLED
#define GET_INDEX(I,J) J,I
#elif   !D1_ENABLED
#error  "Unrecognized interface."
#endif
        ! Open file.
#if     File_ENABLED
#define CLOSE_UNIT close(unit, iostat = err)
        open( file = file & ! LCOV_EXCL_LINE
            , newunit = unit & ! LCOV_EXCL_LINE
            , form = "formatted" & ! LCOV_EXCL_LINE
            , position = "rewind" & ! LCOV_EXCL_LINE
            , access = "sequential" & ! LCOV_EXCL_LINE
            , action = "read" & ! LCOV_EXCL_LINE
            , iostat = err & ! LCOV_EXCL_LINE
            , iomsg = iomsg_def & ! LCOV_EXCL_LINE
            INTEL_SHARED_FILE)
        RETURN_IF_FAILED(__LINE__)
#elif   Unit_ENABLED
#define CLOSE_UNIT
#else
#error  "Unrecognized interface."
#endif
        if (present(roff)) then
            do irow = 1, roff
                read(unit, *, iostat = err, iomsg = iomsg_def)
                RETURN_IF_FAILED(__LINE__)
            end do
        end if
        if (present(header)) then
            call setRecordFrom(unit, header, err, iomsg_def)
            RETURN_IF_FAILED(__LINE__)
        end if
#if     D1_ENABLED && NO_ENABLED
        nrow = RINIT
        call setResized(table, nrow)
        irow = 0_IK
        do
            irow = irow + 1_IK
            if (nrow < irow) then
                nrow = nrow * 2_IK
                call setResized(table, nrow)
            end if
            read(unit, *, iostat = err, iomsg = iomsg_def) table(irow)
            if (err /= 0_IK) then
                if (err == iostat_end) then ! done.
                    if (irow < nrow) call setResized(table, irow - 1_IK)
                    err = 0_IK
                    return
                elseif (present(iomsg)) then
                    iomsg = getStr(__LINE__)//SK_": "//iomsg_def
                end if
                CLOSE_UNIT
                return
            end if
        end do
#elif   D1_ENABLED && TO_ENABLED
        ncol = 0_IK ! We have to determine the number of columns.
        nrow = 1_IK ! We have only one row to read.
        ! Compute the number of table fields.
        blockPresentSep: if (present(sep)) then
            sepLen = len(sep, IK)
            if (sepLen < 1_IK) then
                ncol = 1_IK
                ! This is only either one field or one column.
                exit blockPresentSep
            end if
            if (sep == SK_"," .or. sep == SK_" ") exit blockPresentSep ! .and. sep /= HT ! file can be handled by the Fortran list-directed IO.
            ! \todo The following approach to sep counting must be replaced with a new function like `getFieldSep()` that excludes separators in fields.
            call setRecordFrom(unit, record, err, iomsg_def, 1_IK, lenrec)
            RETURN_IF_FAILED(__LINE__)
            backspace(unit)
            nsep = getCountLoc(record, sep)
            if (nsep == 0_IK) exit blockPresentSep ! can be handled by list-directed IO.
            ncol = nsep + 1_IK
#if         CK_ENABLED
            ncol = ncol / 2_IK
            if (ncol * 2_IK /= nsep + 1_IK) then
                ! the values are not pairs of real and imaginary components.
                if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getErrTableRead(): The number of columns of a complex table must be even."
                err = -1_IK
                CLOSE_UNIT
                return
            end if
            call setResized(field, nsep + 1_IK)
#endif
            call setResized(sepLoc, nsep) ! Pre-allocate the locations of the separators in the record.
            call setResized(table, ncol) ! Initial best guess table size.
            call setRecordFrom(unit, record, err, iomsg_def, 1_IK, lenrec)
            if (err /= 0_IK) then
                if (err == iostat_end) then ! done.
                    err = 0_IK
                elseif (present(iomsg)) then
                    iomsg = getStr(__LINE__)//SK_": "//iomsg_def
                end if
                CLOSE_UNIT
                return
            end if
            call setLoc(sepLoc, lenLoc, record, sep, blindness = sepLen)
            if (lenLoc /= nsep) then
                if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getErrTableRead(): The row "//getStr(irow)// & ! LCOV_EXCL_LINE
                SK_" of the table does not contain the same number of fields as the previous rows."
                err = -1_IK ! The row `irow` of the table rows does not contain `ncol` fields.
                CLOSE_UNIT
                return
            end if
            ! read fields.
#if         CK_ENABLED
            read(record(1 : sepLoc(1) - 1), *, iostat = err, iomsg = iomsg_def) field(1)
            RETURN_IF_FAILED(__LINE__)
            do icol = 2, nsep
                read(record(sepLoc(icol - 1) + sepLen : sepLoc(icol) - 1), *, iostat = err, iomsg = iomsg_def) field(icol)
                RETURN_IF_FAILED(__LINE__)
            end do
            read(record(sepLoc(icol - 1) + sepLen : lenrec), *, iostat = err, iomsg = iomsg_def) field(icol)
            RETURN_IF_FAILED(__LINE__)
            table(1 : ncol) = cmplx(field(1 : nsep : 2), field(2 : nsep + 1 : 2), CKC)
#else
            read(record(1 : sepLoc(1) - 1), *, iostat = err, iomsg = iomsg_def) table(1)
            RETURN_IF_FAILED(__LINE__)
            do icol = 2, nsep
                read(record(sepLoc(icol - 1) + sepLen : sepLoc(icol) - 1), *, iostat = err, iomsg = iomsg_def) table(icol)
                RETURN_IF_FAILED(__LINE__)
            end do
            read(record(sepLoc(icol - 1) + sepLen : lenrec), *, iostat = err, iomsg = iomsg_def) table(icol)
            RETURN_IF_FAILED(__LINE__)
#endif
            return
        end if blockPresentSep
        if (ncol == 0_IK) then
            ! separator can be likely handled by list-directed IO.
#if         SK_ENABLED || CK_ENABLED
            ! Get the separator while respecting quotations.
            record = getFieldSep(unit, SK_", ", fld, ncol, iomsg = iomsg)
#elif       IK_ENABLED || LK_ENABLED || RK_ENABLED
            ! Get the separator.
            record = getFieldSep(unit, SK_", ", ncol, iomsg = iomsg)
#else
#error      "Unrecognized interface."
#endif
        end if
        if (0_IK < ncol) then
            nrow = RINIT
#if         CK_ENABLED
            ! Ensure complex values are parenthesis-delimited.
            call setRecordFrom(unit, record, err, iomsg_def, 1_IK, lenrec)
            if (err /= 0_IK) then
                if (err == iostat_end) then ! done.
                    call setResized(table, 0_IK)
                    err = 0_IK
                elseif (present(iomsg)) then
                    iomsg = getStr(__LINE__)//SK_": "//iomsg_def
                end if
                CLOSE_UNIT
                return
            end if
            backspace(unit)
            irow = getCountLoc(record, SK_"(")
            icol = getCountLoc(record, SK_")")
            if (0_IK == irow .and. 0_IK == icol) then
                ! read the complex table as a simple table of `real` fields.
                !block
                !    real(RKC), allocatable :: rtable(:,:)
                !    err = getErrTableRead(rtable, unit, trans)
                !    return_if_failed
                !    if (present(trans)) then
                !        if (trans) then
                !            do irow = 1, size(rtable,
                !            end do
                !        end if
                !    end if
                !end block
                ! This is not fld format, perhaps csv or similar. Read the table as real in fld format.
                nsep = ncol * 2_IK ! place holder for the number of fields.
                call setResized(field, nsep)
                call setResized(table, ncol) ! Initial best guess table size.
                irow = 1_IK
                read(unit, *, iostat = err, iomsg = iomsg_def) field
                if (err == iostat_end) then
                    err = 0_IK
                    CLOSE_UNIT
                    return
                end if
                RETURN_IF_FAILED(__LINE__)
                table(1 : ncol) = cmplx(field(1 : nsep : 2), field(2 : nsep : 2), CKC)
                return
            elseif (irow /= icol) then
                if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getErrTableRead(): The number of left and right parenthesis delimiters for `complex` table fields must match. '(', ')' = "//getStr([irow, icol])
                err = -1_IK
                CLOSE_UNIT
                return
            end if
            ! The complex table is delimited by `()`. Continue below to read the complex table via Fortran list-directed IO.
#endif
            ! Read the complex table via Fortran list-directed IO.
            call setResized(table, ncol) ! Initial best guess table size.
            irow = 1_IK
            read(unit, *, iostat = err, iomsg = iomsg_def) table(1 : ncol)
            if (err == iostat_end) then
                if (irow < nrow) call setResized(table, ncol)
                err = 0_IK
                CLOSE_UNIT
                return
            end if
            RETURN_IF_FAILED(__LINE__)
        end if
        err = -1_IK
        CLOSE_UNIT
#elif   D2_ENABLED
        ncol = 0_IK ! Assume there is no column in table for now. This is an important assumption.
        ! Compute the number of table fields.
        blockPresentSep: if (present(sep)) then
            sepLen = len(sep, IK)
            if (sepLen < 1_IK) then
                ncol = 1_IK
                ! This is only either one field or one column.
                exit blockPresentSep
            end if
            if (sep == SK_"," .or. sep == SK_" ") exit blockPresentSep ! .and. sep /= HT ! file can be handled by the Fortran list-directed IO.
            ! \todo The following approach to sep counting must be replaced with a new function like `getFieldSep()` that excludes separators in fields.
            call setRecordFrom(unit, record, err, iomsg_def, 1_IK, lenrec)
            RETURN_IF_FAILED(__LINE__)
            backspace(unit)
            nsep = getCountLoc(record, sep)
            if (nsep == 0_IK) exit blockPresentSep
            ncol = nsep + 1_IK
            nrow = RINIT
#if         CK_ENABLED
            ncol = ncol / 2_IK
            if (ncol * 2_IK /= nsep + 1_IK) then
                ! the values are not pairs of real and imaginary components.
                if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getErrTableRead(): The number of columns of a complex table must be even."
                err = -1_IK
                CLOSE_UNIT
                return
            end if
            call setResized(field, nsep + 1_IK)
#endif
            call setResized(sepLoc, nsep) ! Pre-allocate the locations of the separators in the record.
            call setResized(table, [GET_INDEX(nrow, ncol)]) ! Initial best guess table size.
            irow = 0_IK
            loopReadTableRecord: do
                irow = irow + 1_IK
                if (nrow < irow) then
                    nrow = nrow * 2_IK
                    call setResized(table, [GET_INDEX(nrow, ncol)])
                end if
                call setRecordFrom(unit, record, err, iomsg_def, 1_IK, lenrec)
                if (err /= 0_IK) then
                    if (err == iostat_end) then ! done.
                        if (irow < nrow) call setResized(table, [GET_INDEX(irow - 1_IK, ncol)])
                        err = 0_IK
                    elseif (present(iomsg)) then
                        iomsg = getStr(__LINE__)//SK_": "//iomsg_def
                    end if
                    CLOSE_UNIT
                    return
                end if
                call setLoc(sepLoc, lenLoc, record, sep, blindness = sepLen)
                if (lenLoc /= nsep) then
                    if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getErrTableRead(): The row "//getStr(irow)// & ! LCOV_EXCL_LINE
                    SK_" of the table does not contain the same number of fields as the previous rows."
                    err = -1_IK ! The row `irow` of the table rows does not contain `ncol` fields.
                    CLOSE_UNIT
                    return
                end if
                ! read fields.
#if             CK_ENABLED
                read(record(1 : sepLoc(1) - 1), *, iostat = err, iomsg = iomsg_def) field(1)
                RETURN_IF_FAILED(__LINE__)
                do icol = 2, nsep
                    read(record(sepLoc(icol - 1) + sepLen : sepLoc(icol) - 1), *, iostat = err, iomsg = iomsg_def) field(icol)
                    RETURN_IF_FAILED(__LINE__)
                end do
                read(record(sepLoc(icol - 1) + sepLen : lenrec), *, iostat = err, iomsg = iomsg_def) field(icol)
                RETURN_IF_FAILED(__LINE__)
                table(GET_INDEX(irow, 1 : ncol)) = cmplx(field(1 : nsep : 2), field(2 : nsep + 1 : 2), CKC)
#else
                read(record(1 : sepLoc(1) - 1), *, iostat = err, iomsg = iomsg_def) table(GET_INDEX(irow, 1))
                RETURN_IF_FAILED(__LINE__)
                do icol = 2, nsep
                    read(record(sepLoc(icol - 1) + sepLen : sepLoc(icol) - 1), *, iostat = err, iomsg = iomsg_def) table(GET_INDEX(irow, icol))
                    RETURN_IF_FAILED(__LINE__)
                end do
                read(record(sepLoc(icol - 1) + sepLen : lenrec), *, iostat = err, iomsg = iomsg_def) table(GET_INDEX(irow, icol))
                RETURN_IF_FAILED(__LINE__)
#endif
            end do loopReadTableRecord
            ! the flow should never get here.
#if         CHECK_ENABLED
            ! This internal testing can be removed in future.
            error stop MODULE_NAME//SK_"@getErrTableRead(): This is an internal library error."//NLC// & ! LCOV_EXCL_LINE
            SK_"The control flow should never reach this point."//NLC// & ! LCOV_EXCL_LINE
            SK_"Please report this error among with circumstance to the ParaMonte library developers at:"//NLC// & ! LCOV_EXCL_LINE
            SK_"https://github.com/cdslaborg/paramonte/issues."//NLC ! LCOV_EXCL_LINE
#endif
        end if blockPresentSep
        if (ncol == 0_IK) then
            ! separator can be likely handled by list-directed IO.
#if         SK_ENABLED || CK_ENABLED
            ! Get the separator while respecting quotations.
            record = getFieldSep(unit, SK_", ", fld, ncol, iomsg = iomsg)
#elif       IK_ENABLED || LK_ENABLED || RK_ENABLED
            ! Get the separator.
            record = getFieldSep(unit, SK_", ", ncol, iomsg = iomsg)
#else
#error      "Unrecognized interface."
#endif
        end if
        if (0_IK < ncol) then
            nrow = RINIT
#if         CK_ENABLED
            ! Ensure complex values are parenthesis-delimited.
            call setRecordFrom(unit, record, err, iomsg_def, 1_IK, lenrec)
            if (err /= 0_IK) then
                if (err == iostat_end) then ! done.
                    call setResized(table, [0_IK, 0_IK])
                    err = 0_IK
                elseif (present(iomsg)) then
                    iomsg = getStr(__LINE__)//SK_": "//iomsg_def
                end if
                CLOSE_UNIT
                return
            end if
            backspace(unit)
            irow = getCountLoc(record, SK_"(")
            icol = getCountLoc(record, SK_")")
            if (0_IK == irow .and. 0_IK == icol) then
                ! read the complex table as a simple table of `real` fields.
                !block
                !    real(RKC), allocatable :: rtable(:,:)
                !    err = getErrTableRead(rtable, unit, trans)
                !    return_if_failed
                !    if (present(trans)) then
                !        if (trans) then
                !            do irow = 1, size(rtable,
                !            end do
                !        end if
                !    end if
                !end block
                ! This is not fld format, perhaps csv or similar. Read the table as real in fld format.
                nsep = ncol * 2_IK ! place holder for the number of fields.
                call setResized(field, nsep)
                call setResized(table, [GET_INDEX(nrow, ncol)]) ! Initial best guess table size.
                irow = 0_IK
                do
                    irow = irow + 1
                    read(unit, *, iostat = err, iomsg = iomsg_def) field
                    if (err == iostat_end) then
                        if (irow < nrow) call setResized(table, [GET_INDEX(irow - 1_IK, ncol)])
                        err = 0_IK
                        CLOSE_UNIT
                        return
                    end if
                    RETURN_IF_FAILED(__LINE__)
                    table(GET_INDEX(irow, 1 : ncol)) = cmplx(field(1 : nsep : 2), field(2 : nsep : 2), CKC)
                    if (irow < nrow) cycle
                    nrow = nrow * 2_IK
                    call setResized(table, [GET_INDEX(nrow, ncol)])
                end do
            elseif (irow /= icol) then
                if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getErrTableRead(): The number of left and right parenthesis delimiters for `complex` table fields must match. '(', ')' = "//getStr([irow, icol])
                err = -1_IK
                CLOSE_UNIT
                return
            end if
            ! The complex table is delimited by `()`. Continue below to read the complex table via Fortran list-directed IO.
#endif
            ! Read the complex table via Fortran list-directed IO.
            call setResized(table, [GET_INDEX(nrow, ncol)]) ! Initial best guess table size.
            irow = 0_IK
            do
                irow = irow + 1
                read(unit, *, iostat = err, iomsg = iomsg_def) table(GET_INDEX(irow, 1 : ncol))
                if (err == iostat_end) then
                    if (irow < nrow) call setResized(table, [GET_INDEX(irow - 1_IK, ncol)])
                    err = 0_IK
                    CLOSE_UNIT
                    return
                end if
                RETURN_IF_FAILED(__LINE__)
                if (irow < nrow) cycle
                nrow = nrow * 2_IK
                call setResized(table, [GET_INDEX(nrow, ncol)])
            end do
        end if
        err = -1_IK
        CLOSE_UNIT
#endif
#undef  CLOSE_UNIT
#undef  GET_INDEX
#undef  GET_FIELD

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   getErrTableWrite_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     File_ENABLED
        integer(IK) :: unit
#endif
        !character(*, SK), parameter :: gform = SK_"(*(g0,:,','))"
        character(:, SK), allocatable :: format
        character(LEN_IOMSG, SK) :: iomsg_def
#if     D2_ENABLED
        integer(IK) :: nrow, ncol
#endif
        integer(IK) :: irow
#if     File_ENABLED
        if (present(file)) then
            open( file = file & ! LCOV_EXCL_LINE
                , newunit = unit & ! LCOV_EXCL_LINE
                , form = "formatted" & ! LCOV_EXCL_LINE
                , position = "rewind" & ! LCOV_EXCL_LINE
                , access = "sequential" & ! LCOV_EXCL_LINE
                , action = "write" & ! LCOV_EXCL_LINE
                , iostat = err & ! LCOV_EXCL_LINE
                , iomsg = iomsg_def & ! LCOV_EXCL_LINE
                INTEL_SHARED_FILE)
            RETURN_IF_FAILED(__LINE__)
        end if
#elif   Unit_ENABLED
#else
#error  "Unrecognized interface."
#endif
        ! Set the number of columns in format and the quotation formatting. Add `sp,` before NCOL in format to write all numbers in signed format.
#if     NO_ENABLED && D1_ENABLED
#if     CK_ENABLED
#define GET_FORMAT(DELIML,DELIMR,SEP) \
SK_'('//getStrQuoted(DELIML)//SK_',g0,'//SEP//SK_',g0,'//getStrQuoted(DELIMR)//SK_')';
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || RK_ENABLED
#define GET_FORMAT(DELIML,DELIMR,SEP) \
SK_'('//getStrQuoted(DELIML)//SK_',g0,'//getStrQuoted(DELIMR)//SK_')';
#else
#error  "Unrecognized interface."
#endif
#else
#if     CK_ENABLED
#define GET_FORMAT(DELIML,DELIMR,SEP) \
SK_'(*('//getStrQuoted(DELIML)//SK_',g0,'//SEP//SK_',g0,'//getStrQuoted(DELIMR)//SK_',:,'//SEP//SK_'))';
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || RK_ENABLED
#define GET_FORMAT(DELIML,DELIMR,SEP) \
SK_'(*('//getStrQuoted(DELIML)//SK_',g0,'//getStrQuoted(DELIMR)//SK_',:,'//SEP//SK_'))';
#else
#error  "Unrecognized interface."
#endif
#endif
        if (present(sep)) then
            if (present(deliml) .and. present(delimr)) then
                format = GET_FORMAT(deliml,delimr,getStrQuoted(sep))
            elseif (present(deliml)) then
                format = GET_FORMAT(deliml,deliml,getStrQuoted(sep))
            elseif (present(delimr)) then
                format = GET_FORMAT(delimr,delimr,getStrQuoted(sep))
            else
                format = GET_FORMAT(SK_"",SK_"",getStrQuoted(sep))
            end if
        elseif (present(deliml) .and. present(delimr)) then
            format = GET_FORMAT(deliml,delimr,SK_"','")
        elseif (present(deliml)) then
            format = GET_FORMAT(deliml,deliml,SK_"','")
        elseif (present(delimr)) then
            format = GET_FORMAT(delimr,delimr,SK_"','")
        else
            format = GET_FORMAT(SK_"",SK_"",SK_"','")
        end if
        ! Skip lines.
        if (present(roff)) then
            do irow = 1, roff
                write(unit, "(g0)", iostat = err, iomsg = iomsg_def)
                RETURN_IF_FAILED(__LINE__)
            end do
        end if
        ! Define the transposition rules.
#if     NO_ENABLED && D2_ENABLED
        nrow = size(table, 1, IK)
        ncol = size(table, rank(table), IK)
#define TABLE_ROW(I,J) table(I,J)
#elif   TO_ENABLED && D2_ENABLED
        ncol = size(table, 1, IK)
        nrow = size(table, rank(table), IK)
#define TABLE_ROW(I,J) table(J,I)
#elif   !D1_ENABLED
#error  "Unrecognized interface."
#endif
        ! Write header.
        if (present(header)) write(unit, "(g0)", iostat = err, iomsg = iomsg_def) header
        RETURN_IF_FAILED(__LINE__)
        ! Write table.
#if     D1_ENABLED
        write(unit, format, iostat = err, iomsg = iomsg_def) table
        RETURN_IF_FAILED(__LINE__)
#elif   D2_ENABLED
        do irow = 1, nrow
            write(unit, format, iostat = err, iomsg = iomsg_def) TABLE_ROW(irow, 1 : ncol)
            RETURN_IF_FAILED(__LINE__)
        end do
#else
#error  "Unrecognized interface."
#endif
#if     File_ENABLED
        close(unit, iostat = err)
#endif
#undef  GET_FORMAT
#undef  TABLE_ROW

        !%%%%%%%%%%%%%%%%%%
#elif   getFieldSep_ENABLED
        !%%%%%%%%%%%%%%%%%%

        character(LEN_IOMSG, SK) :: iomsg_def
        character(:, SK), allocatable :: record
        integer(IK), parameter :: nsam = 2_IK ! maximum number of line samples.
#if     FFLD_ENABLED
        logical(LK) :: isDiscrete
#endif
#if     File_ENABLED
        integer(IK) :: unit
#elif   !Unit_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define vector `seps` properties.
#if     ID0_ENABLED
#define GET_SIZE(X) len(X, IK)
#define GET_SEPS(I) seps(I:I)
#define OFFSET(I) 0_IK
#elif   CD1_ENABLED
        integer(IK) :: offset(size(seps, 1, IK))
#define GET_SIZE(X) size(X, 1, IK)
#define GET_SEPS(I) seps(I)%val
#define OFFSET(I) offset(I)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: isam, isep, ub, iostat, nseps, freq(GET_SIZE(seps)), freqold(GET_SIZE(seps))
        logical(LK) :: canBeSep(GET_SIZE(seps))
#if     ID0_ENABLED
        nseps = GET_SIZE(seps)
        if (0_IK == nseps) then
            if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getFieldSep(): The condition `0 < len(seps)` must hold."
            sep = SKC_""
            return
        end if
#elif   CD1_ENABLED
        nseps = GET_SIZE(seps)
        if (0_IK < nseps) then
            do isep = 1, nseps
                offset(isep) = len(seps(isep)%val, IK) - 1_IK
                if (offset(isep) < 0_IK) then
                    if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getFieldSep(): The component seps("//getStr(isep)//SK_")%val must have a non-zero length."
                    sep = SKC_""
                    return
                end if
            end do
        else
            if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getFieldSep(): The condition `0 < size(seps)` must hold."
            sep = SKC_""
            return
        end if
#endif
#if     NF_ENABLED
        nfield = 0_IK
#endif
        ! Open file.
#if     File_ENABLED
#define FAILED_RETURN(LINE) \
close(unit, iostat = iostat); if (present(iomsg)) iomsg = getStr(LINE)//SK_": "//iomsg_def; sep = SKC_""; return
        open( file = file & ! LCOV_EXCL_LINE
            , newunit = unit & ! LCOV_EXCL_LINE
            , form = "formatted" & ! LCOV_EXCL_LINE
            , position = "rewind" & ! LCOV_EXCL_LINE
            , access = "sequential" & ! LCOV_EXCL_LINE
            , action = "read" & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
            , iostat = iostat & ! LCOV_EXCL_LINE
            , iomsg = iomsg_def & ! LCOV_EXCL_LINE
            INTEL_SHARED_FILE)
        if (iostat /= 0_IK) then
            if (present(iomsg)) iomsg = iomsg_def
            sep = SKC_""
            return
        end if
#elif   Unit_ENABLED
#define FAILED_RETURN(LINE) \
do isep = 1, isam - 1; backspace(unit); end do; if (present(iomsg)) iomsg = getStr(LINE)//SK_": "//iomsg_def; sep = SKC_""; return ! LCOV_EXCL_LINE
#endif
#if     FDEF_ENABLED
        loopLineSample: do isam = 1, nsam
            freq = 0_IK
            call setRecordFrom(unit, record, iostat, iomsg = iomsg_def, ub = ub)!, linefed = .true._LK)
            !write(*,"(A)") trim(record)
            if (iostat == 0_IK) then
                do isep = 1, nseps
                    freq(isep) = getCountLoc(record(1 : ub), GET_SEPS(isep), blindness = 1_IK + OFFSET(isep))
                end do
            else ! incomplete multiline quote or other reading error.
                FAILED_RETURN(__LINE__)
            end if
            if (1_IK < isam) then
                do isep = 1, nseps
                    canBeSep(isep) = canBeSep(isep) .and. freq(isep) == freqold(isep)
                end do
                freqold = freq
            else
                freqold = freq
                canBeSep = .true._LK
            end if
        end do loopLineSample
#elif   FCSV_ENABLED || FFLD_ENABLED
        block
            character(1,SKC) :: qbeg
            logical(LK) :: quoted
            integer(IK) :: lb, i
            loopLineSample: do isam = 1, nsam
                freq = 0_IK
                qbeg = SKC_" "
                quoted = .false.
                loopReadMultiLineRecord: do
                    lb = 1_IK
                    call setRecordFrom(unit, record, iostat, iomsg_def, lb, ub)!, linefed = .true._LK)
                    if (iostat == 0_IK) then
                        loopSkipQuote: do
                            if (quoted) then
                                ! Skip the quoted field.
                                loopQuoteClose: do i = lb, ub
                                    if (record(i:i) == qbeg) exit loopQuoteClose
                                end do loopQuoteClose
                                quoted = ub < i
                                ! What a nasty quoted field with new line character.
                                if (quoted) cycle loopReadMultiLineRecord ! cycle to read the rest of the field in the next line.
                                i = i + 1_IK
                            else
                                i = lb
                            end if
#if                         FFLD_ENABLED
                            isDiscrete = .true._LK
#endif
                            loopQuoteOpen: do i = i, ub
#if                             FFLD_ENABLED
                                 ! Take care of complex pair first.
                                quoted = record(i:i) == SKC_'('
                                if (quoted) then
                                    lb = i + 1_IK
                                    qbeg = SKC_")"
                                    cycle loopSkipQuote
                                end if
#elif                           !FCSV_ENABLED
#error                          "Unrecognized interface."
#endif
                                quoted = record(i:i) == SKC_"""" .or. record(i:i) == SKC_''''
                                if (quoted) then
                                    lb = i + 1_IK
                                    qbeg = record(i:i)
                                    cycle loopSkipQuote
                                else ! find the separator instances
                                    loopOverSeps: do isep = 1, nseps
                                        if (ub < i + OFFSET(isep)) cycle loopOverSeps
#if                                     FCSV_ENABLED
                                        if (record(i : i + OFFSET(isep)) == GET_SEPS(isep)) freq(isep) = freq(isep) + 1_IK
#elif                                   FFLD_ENABLED
                                        if (record(i : i + OFFSET(isep)) == GET_SEPS(isep)) then
                                            if (isDiscrete .or. GET_SEPS(isep) /= SKC_" ") then
                                                freq(isep) = freq(isep) + 1_IK
                                                isDiscrete = .false._LK
                                            end if
                                        else
                                            isDiscrete = .true._LK
                                        end if
#endif
                                    end do loopOverSeps
                                end if
                            end do loopQuoteOpen
                            exit loopReadMultiLineRecord
                        end do loopSkipQuote
                    end if
                    ! incomplete multiline quote or other reading error.
                    FAILED_RETURN(__LINE__)
                end do loopReadMultiLineRecord
                if (1_IK < isam) then
                    do isep = 1, nseps
                        canBeSep(isep) = canBeSep(isep) .and. freq(isep) == freqold(isep)
                    end do
                    freqold = freq
                else
                    freqold = freq
                    canBeSep = .true._LK
                end if
            end do loopLineSample
        end block
#else
#error  "Unrecognized interface."
#endif
#if     File_ENABLED
        close(unit, iostat = iostat)
        if (iostat /= 0_IK) then
            if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getFieldSep(): Failed to close the input file."
            sep = SKC_""
            return
        end if
#elif   Unit_ENABLED
        ! Backspace the file.
        do isam = 1, nsam
            backspace(unit)
        end do
#endif
        do isep = 1, nseps
            if (canBeSep(isep) .and. 0_IK < freq(isep)) then
#if             NF_ENABLED
                nfield = freq(isep) + 1_IK
#endif
                sep = GET_SEPS(isep)
                return
            end if
        end do
        if (present(iomsg)) iomsg = MODULE_NAME//SK_"@getFieldSep(): There is likely only one field in the records of this file."
        sep = SKC_""
#if     NF_ENABLED
        nfield = 1_IK
#elif   !XX_ENABLED
#error  "Unrecognized interface."
#endif
#undef  OFFSET
#undef  GET_SIZE
#undef  GET_SEPS
#undef  FAILED_RETURN

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getLenFieldMin_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        !use pm_mathNumSys, only: getCountDigit
#if     IK_ENABLED
        ! sign takes one character + `int(log10(huge))` returns
        ! one digit less than the digit count of the huge of actual model.
        lenField = range(mold) + 2_IK
#elif   CK_ENABLED || RK_ENABLED
        !   \bug Intel ifort 2022.
        integer :: rangeLike
        rangeLike = range(mold)
        ! possible leading 0 by some processors, sign, decimal point, exponent symbol, exponent sign
        lenField = precision(mold) + getCountDigit(rangeLike) + 5_IK
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%
#elif   getFormat_ENABLED
        !%%%%%%%%%%%%%%%%

        integer(IK) :: width_def, ndigit_def, subcount_def, lenexp_def, isub
        character(:, SK), allocatable :: field, prefix_def, ed_def, sep_def, deliml_def, subsep_def, delimr_def, sign_str, width_str, count_str, ndigit_str, lenexp_str

        if (getOption(.false._LK, signed)) then
            sign_str = SK_"sp,"
        else
            sign_str = SK_""
        end if

        if (present(prefix)) then
            if (len(prefix) > 0) then
                prefix_def = getStrQuoted(prefix)//SK_"," !//SK_"X,"
            else
                prefix_def = SK_""
            end if
        else
            prefix_def = SK_""
        end if

        if (present(count)) then
            count_str = getStr(count)
            CHECK_ASSERTION(__LINE__, 0_IK < count, SK_": The condition `0 < count` must hold. count = "//getStr(count))
        else
            count_str = SK_"*"
        end if

        if (present(subcount)) then
            subcount_def = subcount
            CHECK_ASSERTION(__LINE__, 0_IK <= subcount, SK_": The condition `0 <= subcount` must hold. subcount = "//getStr(subcount))
        else
#if         CK_ENABLED
            subcount_def = 2_IK
#else
            subcount_def = 1_IK
#endif
        end if

        if (present(ed)) then
            ed_def = getStrLower(ed)
            CHECK_ASSERTION(__LINE__, all([ed_def .in. [character(2,SK) :: 'a', 'e', 'en', 'es', 'ex', 'f', 'g', 'i', 'l']]), SK_": The condition `ed .in. [character(2,SK) :: 'a', 'e', 'en', 'es', 'ex', 'f', 'g', 'i', 'l']` must hold. ed = "//ed)
        else
#if         IK_ENABLED
            ed_def = SK_"i"
#else
            ed_def = SK_"g"
#endif
        end if

        if (present(sep)) then
            sep_def = getStrQuoted(sep)
        else
            sep_def = SK_""", """
        end if

        if (present(subsep)) then
            subsep_def = getStrQuoted(subsep)
        elseif (1_IK < subcount_def) then
            subsep_def = sep_def
        end if

        if (ed_def == SK_"e" .or. ed_def == SK_"en" .or. ed_def == SK_"es" .or. ed_def == SK_"ex" .or. ed_def == SK_"g") then
            if (present(lenexp)) then
                ! per the standard, the exponent field must not be set when g0 is specified or when the precision field is missing.
                CHECK_ASSERTION(__LINE__, 0_IK <= lenexp, SK_": The condition `0 <= lenexp` must hold. lenexp = "//getStr(lenexp))
                !check_assertion(__LINE__, 0_IK < getOption(1_IK, width), SK_": The condition `ed /= 'g' .or. width > 0` must hold.")
                !check_assertion(__LINE__, present(ndigit), SK_": The condition `present(lenexp) .and. present(ndigit) .or. .not. present(lenexp)` must hold.")
                !check_assertion(__LINE__, width_str /= SK_"0", SK_": The condition `present(lenexp) .and. width /= 0 .or. .not. present(lenexp)` must hold.")
                lenexp_str = SK_"e"//getStr(lenexp)
                lenexp_def = lenexp
            else
#if             CK_ENABLED || RK_ENABLED
                lenexp_def = getCountDigit(range(real(0, kind(mold))))
                lenexp_str = SK_"e"//getStr(lenexp_def)
#else
                lenexp_str = SK_""
                lenexp_def = 0_IK
#endif
            end if
        else
            lenexp_str = SK_""
            lenexp_def = 0_IK
        end if

        if (present(ndigit)) then
            CHECK_ASSERTION(__LINE__, 0_IK <= ndigit, SK_": The condition `0 <= ndigit` must hold. ndigit = "//getStr(ndigit))
            ndigit_str = SK_"."//getStr(ndigit)
            ndigit_def = ndigit
        else
#if         CK_ENABLED || RK_ENABLED
            ndigit_def = precision(real(0, kind(mold)))
            ndigit_str = SK_"."//getStr(ndigit_def)
#else
            ndigit_def = 0_IK
            ndigit_str = SK_""
            if (lenexp_str /= SK_"") lenexp_str = SK_""
#endif
        end if

        if (present(width)) then
            width_def = width
            width_str = getStr(width_def)
            ! non-zero width requires non-zero number of digits.
            if (7_IK < width_def .and. ndigit_str == SK_"") then
                ndigit_def = width_def - 7_IK
                ndigit_str = SK_"."//getStr(ndigit_def)
            else
                CHECK_ASSERTION(__LINE__, 0_IK == width .or. (0_IK < width .and. 0_IK < ndigit_def), SK_": The condition `0 == width .or. (0 < width .and. 0 < ndigit)` must hold. width, ndigit = "//getStr([width, ndigit_def]))
            end if
        elseif (ed_def == SK_"a" .or. ed_def == SK_"l") then
            width_def = 0_IK
            width_str = SK_""
#if     !(CK_ENABLED || RK_ENABLED)
        elseif (ed_def == SK_"g") then
            width_def = 0_IK
            width_str = SK_"0"
#endif
        elseif (ed_def == SK_"i") then
            width_def = ndigit_def + 1_IK ! one character for sign.
            width_str = getStr(width_def)
        elseif (ed_def == SK_"f") then
#if         CK_ENABLED || RK_ENABLED
            ! Make it a fixed size field so that it prints nicely on screen.
            !width_def = precision(real(0, kind(mold))) + 3_IK ! three characters for sign, leading 0, and decimal point.
            width_def = precision(real(0, kind(mold))) + ndigit_def + 3_IK ! three characters for sign, leading 0, and decimal point.
            width_str = getStr(width_def)
#else
            width_def = 0_IK ! We do not know the size a priori. Let the compiler set the minimum required size at runtime.
            width_str = getStr(width_def)
#endif
        ! it is a real or complex field with exponent.
        elseif (lenexp_str == SK_"") then ! default field exponent consists of four characters.
            width_def = ndigit_def + 7_IK ! three characters as in `f` descriptor + 4 for default exponent.
            width_str = getStr(width_def)
        elseif (lenexp_def == 0_IK) then ! minimum required field happens when `lenexp = 0` is explicitly specified by the user.
            ! this is tough because we do not know the minimum required exponent length unless we know the type and kind of the field.
            ! therefore, we either use the kind given to us:
#if         CK_ENABLED || RK_ENABLED
            width_def = getCountDigit(range(real(0, kind(mold))))
            width_def = ndigit_def + 5_IK + width_def ! three characters as in `f` descriptor + 2 for exponent symbols + exponent digits.
            width_str = getStr(width_def)
#else
            ! or else, we assume the worst case scenario, the highest precision field of kind \RKB.
            width_def = getCountDigit(range(0._RKB))
            width_def = ndigit_def + 5_IK + width_def ! three characters as in `f` descriptor + 2 for exponent symbols + exponent digits.
            width_str = getStr(width_def)
#endif
        else ! finally, an explicit positive exponent length `lenexp` is known.
            width_def = ndigit_def + 5_IK + lenexp_def ! three characters as in `f` descriptor + 2 for exponent symbols + exponent digits.
            width_str = getStr(width_def)
        end if
        if (width_def == 0_IK .and. lenexp_str /= SK_"") lenexp_str = SK_""

        if (present(deliml) .and. present(delimr)) then
            deliml_def = getStrQuoted(deliml)
            delimr_def = getStrQuoted(delimr)
        elseif (present(deliml)) then
            deliml_def = getStrQuoted(deliml)
            delimr_def = getStrQuoted(deliml)
        elseif (present(delimr)) then
            deliml_def = getStrQuoted(delimr)
            delimr_def = getStrQuoted(delimr)
        else
#if         SK_ENABLED
            deliml_def = getStrQuoted(SK_"""")
            delimr_def = getStrQuoted(SK_"""")
#elif       CK_ENABLED
            deliml_def = getStrQuoted(SK_"(")
            delimr_def = getStrQuoted(SK_")")
#elif       IK_ENABLED || LK_ENABLED || RK_ENABLED || Def_ENABLED
            deliml_def = SK_""
            delimr_def = SK_""
#else
#error      "Unrecognized interface."
#endif
        end if
        if (deliml_def /= SK_"") deliml_def = deliml_def//SK_","
        if (delimr_def /= SK_"") delimr_def = SK_","//delimr_def

        field = ed_def//width_str//ndigit_str//lenexp_str
        if (1_IK < subcount_def) then
            format = field
            do isub = 2, subcount_def
                field = field//SK_","//subsep_def//SK_","//format
            end do
        end if
        format = SK_'('//prefix_def//sign_str//count_str//SK_'('//deliml_def//field//delimr_def//SK_',:,'//sep_def//SK_'))'

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   constructDisplay_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     File_ENABLED
        logical(LK) :: opened
        character(:,SKC), allocatable :: status_def, position_def
        inquire(file = file, opened = opened, number = disp%unit)
        if (opened) close(disp%unit)

        if (present(status)) then
            CHECK_ASSERTION(__LINE__, isValidStatus(status), SK_"@constructDisplay(): The condition `isValidPosition(status)` must hold. status = "//getStr(status))
            status_def = status
        elseif (opened) then
            status_def = SKC_"old"
        else
            status_def = SKC_"unknown"
        end if

        if (present(position)) then
            CHECK_ASSERTION(__LINE__, isValidPosition(position), SK_"@constructDisplay(): The condition `isValidPosition(position)` must hold. position = "//getStr(position))
            position_def = position
        elseif (opened) then
            position_def = SKC_"append"
        else
            position_def = SKC_"asis"
        end if
        !>  \bug
        !>  There is an Intel Fortran compiler bug with the use of `newunit` argument in `open` statement.
        !>  The program opens the file in this procedure. However, it apparently does not keep it open
        !>  in the write methods of the class. Here is the error message:
        !>  forrtl: severe (32): invalid logical unit number, unit -129, file unknown
        !>  update: This could have been due to the finalization routine of the type.
        !disp%unit = getFileUnit()
        open(newunit = disp%unit, file = file, status = status_def, position = position_def)
#elif   Unit_ENABLED
        if (present(unit)) then
            disp%unit = unit
            if (.not. isOpen(unit)) open(unit, status = "scratch")
        else
            disp%unit = output_unit
        end if
#else
#error  "Unrecognized interface."
#endif
        !>  The following setting is critical to prevent closing of the opened file by the `final` subroutine of the class. !disp%opened = .false._LK
        if (present(deliml))    disp%deliml     = deliml
        if (present(delimr))    disp%delimr     = delimr
        if (present(format))    disp%format     = format
        if (present(advance))   disp%advance    = advance
        if (present(tmsize))    disp%tmsize     = tmsize
        if (present(bmsize))    disp%bmsize     = bmsize
        if (present(count))     disp%count      = count

        if (.not. allocated(disp%deliml%string))    disp%deliml%string  = SKC_""
        if (.not. allocated(disp%deliml%integer))   disp%deliml%integer = SKC_""
        if (.not. allocated(disp%deliml%logical))   disp%deliml%logical = SKC_""
        if (.not. allocated(disp%deliml%complex))   disp%deliml%complex = SKC_"("
        if (.not. allocated(disp%deliml%real))      disp%deliml%real    = SKC_""

        if (.not. allocated(disp%delimr%string))    disp%delimr%string  = SKC_""
        if (.not. allocated(disp%delimr%integer))   disp%delimr%integer = SKC_""
        if (.not. allocated(disp%delimr%logical))   disp%delimr%logical = SKC_""
        if (.not. allocated(disp%delimr%complex))   disp%delimr%complex = SKC_")"
        if (.not. allocated(disp%delimr%real))      disp%delimr%real    = SKC_""

        ! Take care of the special cases.
        if (allocated(disp%format%complex)) then
            if (getStrLower(disp%format%complex) == SK_"math") disp%format%complex = FORMAT_GENERIC_DISPLAY_COMPLEX_MATH
        end if

        if (.not. allocated(disp%format%string  ))  disp%format%string  = SKC_'(sp,*('//getStrQuoted(disp%deliml%string  )//SKC_',g0,'           //getStrQuoted(disp%delimr%string  ) //SKC_',:,", "))'
        if (.not. allocated(disp%format%integer ))  disp%format%integer = SKC_'(sp,*('//getStrQuoted(disp%deliml%integer )//SKC_',g0,'           //getStrQuoted(disp%delimr%integer ) //SKC_',:,", "))'
        if (.not. allocated(disp%format%logical ))  disp%format%logical = SKC_'(sp,*('//getStrQuoted(disp%deliml%logical )//SKC_',g0,'           //getStrQuoted(disp%delimr%logical ) //SKC_',:,", "))'
        if (.not. allocated(disp%format%complex ))  disp%format%complex = SKC_'(sp,*('//getStrQuoted(disp%deliml%complex )//SKC_',g0,", ",g0,'   //getStrQuoted(disp%delimr%complex ) //SKC_',:,", "))'
        if (.not. allocated(disp%format%real    ))  disp%format%real    = SKC_'(sp,*('//getStrQuoted(disp%deliml%real    )//SKC_',g0,'           //getStrQuoted(disp%delimr%real    ) //SKC_',:,", "))'

        if (disp%unit == output_unit) then
            if (isFailedGetShellHeight(disp%height)) disp%height = 0_IK
            if (isFailedGetShellWidth(disp%width)) disp%width = 0_IK
        end if

        if (present(text)) then
            disp%text = text
        else
            disp%text = wrap_type(tmsize = disp%tmsize, bmsize = disp%bmsize, width = disp%width, unit = disp%unit, sticky = disp%sticky)
        end if
        if (present(mark)) then
            disp%mark = mark
        else
            disp%mark = mark_type(tmsize = disp%tmsize, bmsize = disp%bmsize, width = disp%width, unit = disp%unit, sticky = disp%sticky)
        end if
        if (present(note)) then
            disp%note = note
        else
            disp%note = note_type(tmsize = disp%tmsize, bmsize = disp%bmsize, width = disp%width, unit = disp%unit, sticky = disp%sticky)
        end if
        if (present(warn)) then
            disp%warn = warn
        else
            disp%warn = warn_type(tmsize = disp%tmsize, bmsize = disp%bmsize, width = disp%width, unit = disp%unit, sticky = disp%sticky)
        end if
        if (present(stop)) then
            disp%stop = stop
        else
            disp%stop = stop_type(tmsize = disp%tmsize, bmsize = disp%bmsize, width = disp%width, unit = disp%unit, sticky = disp%sticky)
        end if
        disp%uninit = .false._LK

        !%%%%%%%%%%%
#elif   show_ENABLED
        !%%%%%%%%%%%

        character(*, SK), parameter     :: NLC = new_line(SK_"a")
        character(:, SK), allocatable   :: def_format
        character(3, SK)                :: def_advance
        integer(IK)                     :: def_tmsize
        integer(IK)                     :: def_bmsize
        integer(IK)                     :: def_count
        integer(IK)                     :: def_unit
        integer(IK)                     :: icount
#if     MEXPRINT_ENABLED
        character(:, SK), allocatable   :: EOL
#endif
        if (self%uninit) then
            select type(self)
            type is (display_type)
                self = display_type()
            end select
        end if
        if (present(sticky)) self%sticky = sticky
        if (present(advance))   then; def_advance   = advance   ; if (self%sticky) self%advance  = advance   ; else; def_advance = self%advance  ; end if
        if (present(tmsize))    then; def_tmsize    = tmsize    ; if (self%sticky) self%tmsize   = tmsize    ; else; def_tmsize  = self%tmsize   ; end if
        if (present(bmsize))    then; def_bmsize    = bmsize    ; if (self%sticky) self%bmsize   = bmsize    ; else; def_bmsize  = self%bmsize   ; end if
        if (present(count))     then; def_count     = count     ; if (self%sticky) self%count    = count     ; else; def_count   = self%count    ; end if
        if (present(unit))      then; def_unit      = unit      ; if (self%sticky) self%unit     = unit      ; else; def_unit    = self%unit     ; end if
#if     MEXPRINT_ENABLED
        call setStrLower(def_advance)
        if (def_advance == "yes") then; EOL = NLC; else; EOL = SK_""; end if
#endif

        ! Set field delimiters and format.
#if     SK_ENABLED
#define FIELD string
#elif   IK_ENABLED
#define FIELD integer
#elif   LK_ENABLED
#define FIELD logical
#elif   CK_ENABLED
#define FIELD complex
#elif   RK_ENABLED
#define FIELD real
#else
#error  "Unrecognized interface."
#endif
#if     CK_ENABLED
#define GET_FORMAT(DELIML,DELIMR) \
SK_'(sp,*('//getStrQuoted(DELIML)//SK_',g0,", ",g0,'//getStrQuoted(DELIMR)//SK_',:,", "))'
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || RK_ENABLED
#define GET_FORMAT(DELIML,DELIMR) \
SK_'(sp,*('//getStrQuoted(DELIML)//SK_',g0,'//getStrQuoted(DELIMR)//SK_',:,", "))'
#else
#error  "Unrecognized interface."
#endif
        if (present(format)) then
            def_format = format
        elseif (present(deliml) .and. present(delimr)) then
            def_format = GET_FORMAT(deliml, delimr)
        elseif (present(deliml)) then
            def_format = GET_FORMAT(deliml, deliml)
        elseif (present(delimr)) then
            def_format = GET_FORMAT(delimr, delimr)
        else
            def_format = self%format%FIELD
        endif
        if (self%sticky) then
            self%format%FIELD = def_format
            if (present(deliml)) self%deliml%FIELD = deliml
            if (present(delimr)) self%delimr%FIELD = delimr
        end if
        ! display contents.
        !   Strategy:
        !   All objects up to rank 2 are directly displayed.
        !   For object os higher rank we recursively reduce the rank by recursively calling the lower rank methods.
        !   The display format follows that of matlab.
        !   That is, by default,
        !       -   A vector is shown as a row.
        !       -   A matrix is shown as a (nrow, ncol) matrix.
        !       -   A cube is shown as a collection of subsequent matrices of shape (:, :, icube).
        !       -   ...
#if     MEXPRINT_ENABLED
#define DISPLAY_NONE \
if (def_unit == output_unit) then; call mexPrintf(NLC); else; write(def_unit, "(A)", advance = def_advance); end if;
#define HEADER(HEAD) \
if (def_unit == output_unit) then; call mexPrintf(HEAD); else; write(def_unit, "(A)", advance = def_advance) HEAD; end if;
#define DISPLAY(ROW) \
if (def_unit == output_unit) then; call mexPrintf(getStr([ROW], def_format)//EOL); else; write(def_unit, def_format, advance = def_advance) ROW; end if
#define MARGIN(SIZE) \
if (def_unit == output_unit) then; call mexPrintf(repeat(NLC, SIZE)); else; write(def_unit, "("//repeat("/", SIZE)//")", advance = "NO"); end if
#else
#define DISPLAY_NONE \
write(def_unit, "(A)", advance = def_advance);
#define HEADER(HEAD) \
write(def_unit, "(A)", advance = def_advance) HEAD;
#define DISPLAY(ROW) \
write(def_unit, def_format, advance = def_advance) ROW;
#define MARGIN(SIZE) \
write(def_unit, "("//repeat("/", SIZE)//")", advance = "NO")
#endif
!#define CALL_DISP(OBJ) call self%show(OBJ, unit = unit, format = def_format, advance = advance, bmsize = 0_IK, tmsize = 0_IK)
        MARGIN(def_tmsize)
        do icount = 1_IK, def_count
#if         CN_ENABLED && D0_ENABLED
            DISPLAY(object)
#elif       CN_ENABLED && D1_ENABLED
            if (0_IK < size(object, 1, IK)) then
                DISPLAY(object)
            else
                DISPLAY_NONE
            end if
#elif       CN_ENABLED && D2_ENABLED
            block
                integer(IK) :: irow
                do irow = 1, size(object, 1, IK)
                    if (0_IK < size(object, 2, IK)) then
                        DISPLAY(object(irow, :))
                    else
                        DISPLAY_NONE
                    end if
                end do
            end block
#elif       CN_ENABLED && D3_ENABLED
            block
                integer(IK) :: imat, irow
                do imat = 1, size(object, 3, IK)
                    HEADER(SK_"slice(:,:,"//getStr(imat)//SK_") = ")
                    do irow = 1, size(object, 1, IK)
                        if (0_IK < size(object, 2, IK)) then
                            DISPLAY(object(irow, :, imat))
                        else
                            DISPLAY_NONE
                        end if
                    end do
                end do
            end block
#elif       (BS_ENABLED || PS_ENABLED) && D0_ENABLED
            DISPLAY(object%val)
#elif       (BS_ENABLED || PS_ENABLED) && D1_ENABLED
            block
                integer(IK) :: idim
                if (0_IK < size(object, 1, IK)) then
                    DISPLAY((object(idim)%val, idim = 1, size(object, 1, IK)))
                else
                    DISPLAY_NONE
                end if
            end block
#elif       (BS_ENABLED || PS_ENABLED) && D2_ENABLED
            block
                integer(IK) :: icol, irow
                do irow = 1, size(object, 1, IK)
                    if (0_IK < size(object, 2, IK)) then
                        DISPLAY((object(irow, icol)%val, icol = 1, size(object, 2, IK)))
                    else
                        DISPLAY_NONE
                    end if
                end do
            end block
#elif       (BS_ENABLED || PS_ENABLED) && D3_ENABLED
            block
                integer(IK) :: imat, irow, icol
                do imat = 1, size(object, 3, IK)
                    HEADER(SK_"object(:,:,"//getStr(imat)//SK_")%val = ")
                    do irow = 1, size(object, 1, IK)
                        if (0_IK < size(object, 2, IK)) then
                            DISPLAY((object(irow, icol, imat)%val, icol = 1, size(object, 2, IK)))
                        else
                            DISPLAY_NONE
                        end if
                    end do
                end do
            end block
#elif       (BV_ENABLED || PV_ENABLED) && D0_ENABLED
            DISPLAY(object%val)
#elif       (BV_ENABLED || PV_ENABLED) && D1_ENABLED
            block
                integer(IK) :: irow
                do irow = 1, size(object, 1, IK)
                    if (0_IK < size(object(irow)%val, 1, IK)) then
                        DISPLAY(object(irow)%val(:))
                    else
                        DISPLAY_NONE
                    end if
                end do
            end block
#elif       (BV_ENABLED || PV_ENABLED) && D2_ENABLED
            block
                integer(IK) :: irow, imat
                do imat = 1, size(object, 2, IK)
                    HEADER(SK_"object(:,"//getStr(imat)//SK_")%val(:) = ")
                    do irow = 1, size(object, 1, IK)
                        if (0_IK < size(object(irow, imat)%val, 1, IK)) then
                            DISPLAY(object(irow, imat)%val(:))
                        else
                            DISPLAY_NONE
                        end if
                    end do
                end do
            end block
#elif       (BM_ENABLED || CM_ENABLED) && D0_ENABLED
            block
                integer(IK) :: irow
                do irow = 1, size(object%val, 1, IK)
                    if (0_IK < size(object%val, 2, IK)) then
                        DISPLAY(object%val(irow, :))
                    else
                        DISPLAY_NONE
                    end if
                end do
            end block
#elif       (BM_ENABLED || CM_ENABLED) && D1_ENABLED
            block
                integer(IK) :: imat, irow
                do imat = 1, size(object, 1, IK)
                    HEADER(SK_"object("//getStr(imat)//SK_")%val(:,:) = ")
                    do irow = 1, size(object(imat)%val, 1, IK)
                        if (0_IK < size(object(imat)%val, 2, IK)) then
                            DISPLAY(object(imat)%val(irow, :))
                        else
                            DISPLAY_NONE
                        end if
                    end do
                end do
            end block
#elif       (BC_ENABLED || PC_ENABLED) && D0_ENABLED
            block
                integer(IK) :: imat, irow
                do imat = 1, size(object%val, 3, IK)
                    HEADER(SK_"object%val(:,:,"//getStr(imat)//SK_") = ")
                    do irow = 1, size(object%val, 1, IK)
                        if (0_IK < size(object%val, 2, IK)) then
                            DISPLAY(object%val(irow, :, imat))
                        else
                            DISPLAY_NONE
                        end if
                    end do
                end do
            end block
#else
#error      "Unrecognized interface."
#endif
        end do
        MARGIN(def_bmsize)
        flush(def_unit)
#undef  DISPLAY_NONE
#undef  GET_FORMAT
#undef  DISPLAY
#undef  FIELD

        !%%%%%%%%%%%
#elif   dump_ENABLED
        !%%%%%%%%%%%

#define CALL_SHOW \
call self%show(object, tmsize = tmsize, bmsize = bmsize, count = count, unit = unit, format = format, advance = advance, sticky = sticky)
        select type (object)
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if         SK5_ENABLED
            type is (character(*,SK5)); CALL_SHOW
#endif
#if         SK4_ENABLED
            type is (character(*,SK4)); CALL_SHOW
#endif
#if         SK3_ENABLED
            type is (character(*,SK3)); CALL_SHOW
#endif
#if         SK2_ENABLED
            type is (character(*,SK2)); CALL_SHOW
#endif
#if         SK1_ENABLED
            type is (character(*,SK1)); CALL_SHOW
#endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if         IK5_ENABLED
            type is (integer(IK5)); CALL_SHOW
#endif
#if         IK4_ENABLED
            type is (integer(IK4)); CALL_SHOW
#endif
#if         IK3_ENABLED
            type is (integer(IK3)); CALL_SHOW
#endif
#if         IK2_ENABLED
            type is (integer(IK2)); CALL_SHOW
#endif
#if         IK1_ENABLED
            type is (integer(IK1)); CALL_SHOW
#endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if         LK5_ENABLED
            type is (logical(LK5)); CALL_SHOW
#endif
#if         LK4_ENABLED
            type is (logical(LK4)); CALL_SHOW
#endif
#if         LK3_ENABLED
            type is (logical(LK3)); CALL_SHOW
#endif
#if         LK2_ENABLED
            type is (logical(LK2)); CALL_SHOW
#endif
#if         LK1_ENABLED
            type is (logical(LK1)); CALL_SHOW
#endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if         CK5_ENABLED
            type is (complex(CK5)); CALL_SHOW
#endif
#if         CK4_ENABLED
            type is (complex(CK4)); CALL_SHOW
#endif
#if         CK3_ENABLED
            type is (complex(CK3)); CALL_SHOW
#endif
#if         CK2_ENABLED
            type is (complex(CK2)); CALL_SHOW
#endif
#if         CK1_ENABLED
            type is (complex(CK1)); CALL_SHOW
#endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if         RK5_ENABLED
            type is (real(RK5)); CALL_SHOW
#endif
#if         RK4_ENABLED
            type is (real(RK4)); CALL_SHOW
#endif
#if         RK3_ENABLED
            type is (real(RK3)); CALL_SHOW
#endif
#if         RK2_ENABLED
            type is (real(RK2)); CALL_SHOW
#endif
#if         RK1_ENABLED
            type is (real(RK1)); CALL_SHOW
#endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            class default; error stop "Unrecognized unsupported type for dumping to display." ! LCOV_EXCL_LINE
        end select
#undef  CALL_SHOW

        !%%%%%%%%%%%
#elif   wrap_ENABLED
        !%%%%%%%%%%%

!#if     __GFORTRAN__
!        ! Bypass gfortran bug with deallocation of heap arrays.
!#define BYPASS_GFORTRAN_BUG if (allocated(temp)) deallocate(temp);
!#else
!#define BYPASS_GFORTRAN_BUG
!#endif
!
!#define RESIZE(strWrapped) \
!        if (pos + fullwidth + lenLF > lenStrWrapped) then; \
!            allocate(character(lenStrWrapped * 2, SKC) :: temp); \
!            temp(1:lenStrWrapped) = strWrapped; \
!            call move_alloc(temp, strWrapped); \
!            lenStrWrapped = len(strWrapped, IK); \
!            BYPASS_GFORTRAN_BUG \
!        end if
       !character(:,SKC), allocatable   :: temp
        character(1,SKC), parameter     :: FILL_SKC = MFILL
        character(1,SKC), parameter     :: LF = new_line(SKC_"a") ! char(10, SKC)
        integer(IK)     , parameter     :: lenLF = len(LF, IK)

        integer(IK)                     :: def_lwsize, def_rwsize, def_twsize, def_bwsize, def_tmsize, def_bmsize, def_width, def_unit
        integer(IK)                     :: lenStr, lenLine, lenNewLine, fullwidth, pos, i, iend !, lenStrWrapped
        character(1,SKC)                :: def_lwfill, def_rwfill, def_twfill, def_bwfill, def_fill
        character(:,SKC), allocatable   :: tbwrap, def_newline, strWrapped

        if (self%uninit) then
            if (.not. allocated(self%newline)) self%newline = LF
            if (.not. allocated(self%lwfill )) self%lwfill  = FILL_SKC
            if (.not. allocated(self%rwfill )) self%rwfill  = FILL_SKC
            if (.not. allocated(self%twfill )) self%twfill  = FILL_SKC
            if (.not. allocated(self%bwfill )) self%bwfill  = FILL_SKC
            if (.not. allocated(self%fill   )) self%fill    = SKC_" "
            self%uninit = .true._LK
        end if

        if (present(sticky)) self%sticky = sticky
        if (present(width   )) then; def_width  = width     ; if (self%sticky) self%width   = width     ; else; def_width   = self%width    ; end if
        if (present(lwfill  )) then; def_lwfill = lwfill    ; if (self%sticky) self%lwfill  = lwfill    ; else; def_lwfill  = self%lwfill   ; end if
        if (present(rwfill  )) then; def_rwfill = rwfill    ; if (self%sticky) self%rwfill  = rwfill    ; else; def_rwfill  = self%rwfill   ; end if
        if (present(twfill  )) then; def_twfill = twfill    ; if (self%sticky) self%twfill  = twfill    ; else; def_twfill  = self%twfill   ; end if
        if (present(bwfill  )) then; def_bwfill = bwfill    ; if (self%sticky) self%bwfill  = bwfill    ; else; def_bwfill  = self%bwfill   ; end if
        if (present(lwsize  )) then; def_lwsize = lwsize    ; if (self%sticky) self%lwsize  = lwsize    ; else; def_lwsize  = self%lwsize   ; end if
        if (present(rwsize  )) then; def_rwsize = rwsize    ; if (self%sticky) self%rwsize  = rwsize    ; else; def_rwsize  = self%rwsize   ; end if
        if (present(twsize  )) then; def_twsize = twsize    ; if (self%sticky) self%twsize  = twsize    ; else; def_twsize  = self%twsize   ; end if
        if (present(bwsize  )) then; def_bwsize = bwsize    ; if (self%sticky) self%bwsize  = bwsize    ; else; def_bwsize  = self%bwsize   ; end if
        if (present(tmsize  )) then; def_tmsize = tmsize    ; if (self%sticky) self%tmsize  = tmsize    ; else; def_tmsize  = self%tmsize   ; end if
        if (present(bmsize  )) then; def_bmsize = bmsize    ; if (self%sticky) self%bmsize  = bmsize    ; else; def_bmsize  = self%bmsize   ; end if
        if (present(fill    )) then; def_fill   = fill      ; if (self%sticky) self%fill    = fill      ; else; def_fill    = self%fill     ; end if
        if (present(unit    )) then; def_unit   = unit      ; if (self%sticky) self%unit    = unit      ; else; def_unit    = self%unit     ; end if
        if (present(newline )) then; def_newline= newline   ; if (self%sticky) self%newline = newline   ; else; def_newline = self%newline  ; end if

        lenStr = len(str, IK)
        lenNewLine = len(def_newline, IK)
        lenLine = lenStr * 5_IK / def_width ! best guess for the number of lines, assuming each line is 1/5 of the specified width.
        fullwidth = def_lwsize + def_width + def_rwsize ! the full width of each line (excluding the newline character at the end).
        !lenStrWrapped = (fullwidth + lenLF) * (def_twsize + def_bwsize + lenLine)
        allocate(character((fullwidth + lenLF) * (def_twsize + def_bwsize + lenLine) + (def_tmsize + def_bmsize) * lenLF, SKC) :: strWrapped)

        ! Add the top margin.

        pos = def_tmsize * lenLF
        strWrapped(1 : pos) = repeat(LF, def_tmsize)
        tbwrap = repeat(def_twfill, fullwidth)//LF
        do i = 1_IK, def_twsize
            strWrapped(pos + 1 : pos + fullwidth + lenLF) = tbwrap
            pos = pos + fullwidth + lenLF
        end do

        i = 1_IK
        iend = 0_IK
        do
            if (iend == lenStr) exit
            lenLine = index(str(i : lenStr), def_newline)
            if (lenLine > 0_IK) then
                iend = i + lenLine - 2_IK
            else
                iend = lenStr
            end if
            !RESIZE(strWrapped)
            if (len(strWrapped, IK) < pos + fullwidth + lenLF) call setResized(strWrapped)
            !write(*, *) "len(str(i:iend))", len(str(i:iend))
            call setCentered( strWrapped(pos + 1_IK : pos + fullwidth) & ! LCOV_EXCL_LINE
                            , str(i : iend) & ! LCOV_EXCL_LINE
                            , lmsize = def_lwsize   & ! LCOV_EXCL_LINE
                            , rmsize = def_rwsize   & ! LCOV_EXCL_LINE
                            , lmfill = def_lwfill   & ! LCOV_EXCL_LINE
                            , rmfill = def_rwfill   & ! LCOV_EXCL_LINE
                            , fill = def_fill       & ! LCOV_EXCL_LINE
                            )
            strWrapped(pos + fullwidth + 1_IK : pos + fullwidth + lenLF) = LF
            pos = pos + fullwidth + lenLF
            i = iend + lenNewLine + 1_IK
        end do

        ! Add the bottom margin.

        tbwrap(1 : fullwidth) = repeat(def_bwfill, fullwidth)

        do i = 1_IK, def_bwsize
            !RESIZE(strWrapped)
            if (pos + fullwidth + lenLF > len(strWrapped, IK)) call setResized(strWrapped)
            strWrapped(pos + 1 : pos + fullwidth) = tbwrap
            pos = pos + fullwidth
        end do

        i = pos + def_bmsize * lenLF
        if (len(strWrapped, IK) < i) call setResized(strWrapped)
        strWrapped(pos + 1 : i) = repeat(LF, def_bmsize)
#if     MEXPRINT_ENABLED
        if (def_unit == output_unit) then
            call mexPrintf(strWrapped(1 : i)//LF)
        else
            write(def_unit, "(a)") strWrapped(1 : i)
        end if
#else
        write(def_unit, "(a)") strWrapped(1 : i)
#endif
        flush(def_unit)
        !write(*, "(a)")
        !write(*, "(a)") strWrapped(1:pos)
        !write(*, "(a)")
        deallocate(strWrapped)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CATCH_ERR_IF_FAILED
#undef  INTEL_SHARED_FILE
#undef  RETURN_IF_FAILED
#undef  IOSTAT_IOMSG
#undef  SET_STAT_IO
#undef  ITEM