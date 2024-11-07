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
!>  This file contains procedure implementations of [pm_io](@ref pm_io).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if MEXPRINT_ENABLED
#include "fintrf.h"
#endif

submodule (pm_io) routines

#if CHECK_ENABLED
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_str, only: NLC
    use pm_kind, only: RKB
    use pm_err, only: getFine
    use pm_err, only: getLine
    use pm_val2str, only: getStr
    use pm_option, only: getOption
    use pm_arrayFind, only: setLoc
    use pm_arrayResize, only: setResized
    use pm_arrayCenter, only: setCentered
    use pm_mathNumSys, only: getCountDigit
    use pm_arrayMembership, only: operator(.in.)
    use pm_sysShell, only: isFailedGetShellHeight
    use pm_sysShell, only: isFailedGetShellWidth
    use iso_fortran_env, only: output_unit
    use pm_mathNumSys, only: getCountDigit
    use pm_strASCII, only: getStrQuoted
    use pm_arrayFind, only: getCountLoc
    use pm_strASCII, only: setStrLower
    use pm_strASCII, only: getStrLower

    use iso_fortran_env, only: iostat_end, iostat_eor

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure openArg_typer

        !character(:, SK), allocatable :: msg
        !if (allocated(msg)) deallocate(msg)

        construction_block: block

            ! Set the file access mode.

            if (present(access)) then
                openArg%access = access
                !if (isFalseAssertion(isValidAccess(openArg%access), SK_"The input `access` value is invalid: "//getStr(openArg%access), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidAccess(openArg%access)) then
                    openArg%iomsg = SK_"The input `access` value is invalid: "//getStr(openArg%access)
                    exit construction_block
                end if
            end if

            ! Set the file action mode.

            if (present(action)) then
                openArg%action = action
                !if (isFalseAssertion(isValidAction(openArg%action), SK_"The input `action` value is invalid: "//getStr(openArg%action), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidAction(openArg%action)) then
                    openArg%iomsg = SK_"The input `action` value is invalid: "//getStr(openArg%action)
                    exit construction_block
                end if
            end if

            ! Set the file asynchronous mode.

            if (present(asynchronous)) then
                openArg%asynchronous = asynchronous
                !if (isFalseAssertion(isValidAsynchronous(openArg%asynchronous), SK_"The input `asynchronous` value is invalid: "//getStr(openArg%asynchronous), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidAsynchronous(openArg%asynchronous)) then
                    openArg%iomsg = SK_"The input `asynchronous` value is invalid: "//getStr(openArg%asynchronous)
                    exit construction_block
                end if
            end if

            ! Set the file number padding style (with blanks or zeros).

            if (present(blank)) then
                openArg%blank = blank
                !if (isFalseAssertion(isValidBlank(openArg%blank), SK_"The input `blank` value is invalid: "//getStr(openArg%blank), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidBlank(openArg%blank)) then
                    openArg%iomsg = SK_"The input `blank` value is invalid: "//getStr(openArg%blank)
                    exit construction_block
                end if
            end if

            ! Set the file decimal symbol.

            if (present(decimal)) then
                openArg%decimal = decimal
                !if (isFalseAssertion(isValidDecimal(openArg%decimal), SK_"The input `decimal` value is invalid: "//getStr(openArg%decimal), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidDecimal(openArg%decimal)) then
                    openArg%iomsg = SK_"The input `decimal` value is invalid: "//getStr(openArg%decimal)
                    exit construction_block
                end if
            end if

            ! Set the file string delimiter (string quotation).

            if (present(delim)) then
                openArg%delim = delim
                !if (isFalseAssertion(isValidDelim(openArg%delim), SK_"The input `delim` value is invalid: "//getStr(openArg%delim), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidDelim(openArg%delim)) then
                    openArg%iomsg = SK_"The input `delim` value is invalid: "//getStr(openArg%delim)
                    exit construction_block
                end if
            end if

            ! Set the file encoding.

            if (present(encoding)) then
                openArg%encoding = encoding
                !if (isFalseAssertion(isValidEncoding(openArg%encoding), SK_"The input `encoding` value is invalid: "//getStr(openArg%encoding), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidEncoding(openArg%encoding)) then
                    openArg%iomsg = SK_"The input `encoding` value is invalid: "//getStr(openArg%encoding)
                    exit construction_block
                end if
            end if

            ! Set the file form.

            if (present(form)) then
                openArg%form = form
                !if (isFalseAssertion(isValidForm(openArg%form), SK_"The input `form` value is invalid: "//getStr(openArg%form), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidForm(openArg%form)) then
                    openArg%iomsg = SK_"The input `form` value is invalid: "//getStr(openArg%form)
                    exit construction_block
                end if
            elseif (openArg%access == SK_"direct" .or. openArg%access == SK_"stream") then
                openArg%form = SK_"unformatted"
            end if

            ! Set the file variable padding mode.

            if (present(pad)) then
                openArg%pad = pad
                !if (isFalseAssertion(isValidPad(openArg%pad), SK_"The input `pad` value is invalid: "//getStr(openArg%pad), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidPad(openArg%pad)) then
                    openArg%iomsg = SK_"The input `pad` value is invalid: "//getStr(openArg%pad)
                    exit construction_block
                end if
            end if

            ! Set the file variable padding mode.

            if (present(position)) then
                openArg%position = position
                !if (isFalseAssertion(isValidPosition(openArg%position), SK_"The input `position` value is invalid: "//getStr(openArg%position), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidPosition(openArg%position)) then
                    openArg%iomsg = SK_"The input `position` value is invalid: "//getStr(openArg%position)
                    exit construction_block
                end if
            end if

            ! Set the file record length.

            if (openArg%access == SK_"sequential") then
                if (present(recl)) then
                    openArg%recl = recl
                    !if (isFalseAssertion(isValidRecl(openArg%recl), SK_"The input `recl` value is invalid: "//getStr(openArg%recl), iostat, iomsg)) return ! LCOV_EXCL_LINE
                    if(.not. isValidRecl(openArg%recl)) then
                        openArg%iomsg = SK_"The input `recl` value is invalid: "//getStr(openArg%recl)
                        exit construction_block
                    end if
                end if
            else
                !if (isFalseAssertion(logical(openArg%access == SK_"stream" .and. .not. present(recl), LK), SK_"The input argument `recl` must not be present when the file access mode is `stream`.", iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. (logical(openArg%access == SK_"stream" .and. .not. present(recl), LK))) then
                    openArg%iomsg = SK_"The input argument `recl` must not be present when the file access mode is `stream`."
                    exit construction_block
                end if
                !if (isFalseAssertion(logical(openArg%access == SK_"direct" .and. .not. present(recl), LK), SK_"The input argument `recl` must be present when the file access mode is `direct`.", iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. (logical(openArg%access == SK_"direct" .and. .not. present(recl), LK))) then
                    openArg%iomsg = SK_"The input argument `recl` must be present when the file access mode is `direct`."
                    exit construction_block
                end if
            end if

            ! Set the file rounding mode.

            if (present(round)) then
                openArg%round = round
                !if (isFalseAssertion(isValidRound(openArg%round), SK_"The input `round` value is invalid: "//getStr(openArg%round), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidRound(openArg%round)) then
                    openArg%iomsg = SK_"The input `round` value is invalid: "//getStr(openArg%round)
                    exit construction_block
                end if
            end if

            ! Set the file sign mode.

            if (present(sign)) then
                openArg%sign = sign
                !if (isFalseAssertion(isValidSign(openArg%sign), SK_"The input `sign` value is invalid: "//getStr(openArg%sign), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidSign(openArg%sign)) then
                    openArg%iomsg = SK_"The input `sign` value is invalid: "//getStr(openArg%sign)
                    exit construction_block
                end if
            end if

            ! Set the file status.

            if (present(status)) then
                openArg%status = status
                !if (isFalseAssertion(isValidStatus(openArg%status), SK_"The input `status` value is invalid: "//getStr(openArg%status), iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(.not. isValidStatus(openArg%status)) then
                    openArg%iomsg = SK_"The input `status` value is invalid: "//getStr(openArg%status)
                    exit construction_block
                end if
            end if

            ! Set the file name.

            if (present(file)) then
                openArg%file = file
                !if (isFalseAssertion(logical(openArg%status == SK_"scratch", LK), SK_"A file name cannot be specified when `status = ""scratch""`", iostat, iomsg)) return ! LCOV_EXCL_LINE
                if(openArg%status == SK_"scratch") then
                    openArg%iomsg = SK_"A file name cannot be specified when `status = ""scratch""`"
                    exit construction_block
                end if
            end if

            ! Set the file unit.

            openArg%unit = getOption(getFileUnit(file), unit)
            if (openArg%unit == -1_IK) openArg%unit = getFileUnit()

            if (present(iostat)) iostat = 0_IK
            return ! this is critical for safe return.

            ! Set the iomsg component length.

            if (present(iomsg)) then
                openArg%iomsg = repeat(" ", len(iomsg, IK))
            else
                openArg%iomsg = repeat(" ", LEN_IOMSG)
            end if

            if (present(iostat)) iostat = 0_IK
            return ! this is critical for safe return.

        end block construction_block

        openArg%iomsg = SK_"FATAL RUNTIME ERROR: "//openArg%iomsg
        if (present(iomsg)) iomsg = openArg%iomsg
        if (present(iostat)) then
            iostat = 1_IK
        else
            error stop openArg%iomsg ! LCOV_EXCL_LINE
        end if

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
    !>  \brief
    !>  Generate and return `.false.` if the input assertion is `.true.`, otherwise
    !>  return `.false.`, all the while setting the rest of the `optional` output arguments.
    !>
    !>  \details
    !>  If the `optional` output arguments and missing and the `assertion` is `.false.`,
    !>  then program will halt by calling the `error stop msg` statement, where `msg` is the input message to display.
    !>
    !>  \param[in]      assertion   :   The input scalar `logical` of default kind \LK representing the logical value of assertion to evaluate.
    !>  \param[in]      msg         :   The input scalar `character` of default kind \SK representing the contents of the output argument `iomsg`.
    !>                                  If the output argument `iostat` is missing, then `msg` will be passed to the `error stop` statement.
    !>  \param[in]      iostat      :   The output scalar `integer` of default kind \IK that is set to `1` if assertion is `.false.`, otherwise `0`.<br>
    !>                                  (**optional**. If missing and `assertion` is .false., then program will halt by calling `error stop`).
    !>  \param[inout]   iomsg       :   The input/output scalar `character` of default kind \SK that contains the input `msg`.<br>
    !>                                  (**optional**. The presence of the output `iomsg` is relevant only if `iostat` is also present.)
    !>
    !>  \return
    !>  `assertionIsFalse`          :   The output scalar `logical` of default kind \LK that is `.true.` if `assertion == .false.` holds, otherwise, it is `.false.`.
    !>
    !>  \interface{isFalseAssertion}
    !>  \code{.F90}
    !>
    !>      use pm_io, only: isFalseAssertion
    !>
    !>      assertionIsFalse = isFalseAssertion(assertion, msg, iostat = iostat, iomsg = iomsg)
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This procedure is primarily intended for handling of IO errors within the parent module of the procedure ([pm_io](@ref pm_io)).
    !>
    !>  \see
    !>  [openArg_type](@ref pm_io::openArg_type)<br>
    !>
    !>  \test
    !>  [test_pm_io](@ref test_pm_io)
    !>
    !>  \final{isFalseAssertion}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 2:22 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function isFalseAssertion(assertion, msg, iostat, iomsg) result(assertionIsFalse)
        #if __INTEL_COMPILER
                !DEC$ ATTRIBUTES INLINE :: isFalseAssertion
        #endif
        logical(LK)     , intent(in)                :: assertion
        character(*, SK), intent(in)                :: msg
        integer(IK)     , intent(out)   , optional  :: iostat
        character(*, SK), intent(inout) , optional  :: iomsg
        logical(LK)                                 :: assertionIsFalse
        assertionIsFalse = .not. assertion
        if (assertionIsFalse .and. present(iostat)) then
            if (present(iostat)) then
                if (present(iomsg)) iomsg = trim(adjustl(msg))
                iostat = 1_IK
            else
                error stop SK_"FATAL RUNTIME ERROR: "//trim(adjustl(msg)) ! LCOV_EXCL_LINE
            end if
        elseif (present(iostat)) then
            iostat = 0_IK
        end if
    end function
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define field_typer_ENABLED 1

    module procedure field_typer
#include "pm_io@routines.inc.F90"
    end procedure

#undef field_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getAction_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    module procedure getActionFile
#include "pm_io@routines.inc.F90"
    end procedure

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    module procedure getActionUnit
#include "pm_io@routines.inc.F90"
    end procedure

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getAction_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isOpen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    module procedure isOpenFile
#include "pm_io@routines.inc.F90"
    end procedure

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    module procedure isOpenUnit
#include "pm_io@routines.inc.F90"
    end procedure

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isOpen_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCountRecord_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    module procedure getCountRecordFile
#include "pm_io@routines.inc.F90"
    end procedure

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    module procedure getCountRecordUnit
#include "pm_io@routines.inc.F90"
    end procedure

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCountRecord_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getContentsFrom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    module procedure getContentsFromFile_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    module procedure getContentsFromUnit_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getContentsFrom_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setContentsFrom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CDD_ENABLED 1

    module procedure setContentsFromFileCDD_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CII_ENABLED 1

    module procedure setContentsFromFileCII_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CDD_ENABLED 1

    module procedure setContentsFromUnitCDD_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CII_ENABLED 1

    module procedure setContentsFromUnitCII_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setContentsFrom_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setContentsTo_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CDD_ENABLED 1

    module procedure setContentsToFileCDD_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CII_ENABLED 1

    module procedure setContentsToFileCII_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define Unit_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CDD_ENABLED 1
!
!    module procedure setContentsToUnitCDD_SK
!        use pm_kind, only: SKG => SK
!#include "pm_io@routines.inc.F90"
!    end procedure
!
!#undef CDD_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CII_ENABLED 1
!
!    module procedure setContentsToUnitCII_SK
!        use pm_kind, only: SKG => SK
!#include "pm_io@routines.inc.F90"
!    end procedure
!
!#undef CII_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef Unit_ENABLED
!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setContentsTo_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setFileClosed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure setFileClosed_IK
        use pm_kind, only: IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setFileClosed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCountRecordLeft_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure getCountRecordLeft_IK
        use pm_kind, only: IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCountRecordLeft_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRecordFrom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    module procedure getRecordFromUnit_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRecordFrom_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRecordFrom_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UR_ENABLED 1

    module procedure setRecordFromUR_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef UR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define URII_ENABLED 1

    module procedure setRecordFromURII_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef URII_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRecordFrom_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getErrTableRead_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadFile_NO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadFile_NO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadFile_NO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadFile_NO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadFile_NO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadFile_NO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadFile_NO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadFile_NO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadFile_NO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadFile_NO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadFile_NO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadFile_NO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadFile_NO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadFile_NO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadFile_NO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadFile_NO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadFile_NO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadFile_NO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadFile_NO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadFile_NO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadFile_NO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadFile_NO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadFile_NO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadFile_NO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadFile_NO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadUnit_NO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadUnit_NO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadUnit_NO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadUnit_NO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadUnit_NO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadUnit_NO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadUnit_NO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadUnit_NO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadUnit_NO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadUnit_NO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadUnit_NO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadUnit_NO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadUnit_NO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadUnit_NO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadUnit_NO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadUnit_NO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadUnit_NO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadUnit_NO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadUnit_NO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadUnit_NO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadUnit_NO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadUnit_NO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadUnit_NO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadUnit_NO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadUnit_NO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadFile_NO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadFile_NO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadFile_NO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadFile_NO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadFile_NO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadFile_NO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadFile_NO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadFile_NO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadFile_NO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadFile_NO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadFile_NO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadFile_NO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadFile_NO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadFile_NO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadFile_NO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadFile_NO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadFile_NO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadFile_NO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadFile_NO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadFile_NO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadFile_NO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadFile_NO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadFile_NO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadFile_NO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadFile_NO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadUnit_NO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadUnit_NO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadUnit_NO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadUnit_NO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadUnit_NO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadUnit_NO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadUnit_NO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadUnit_NO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadUnit_NO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadUnit_NO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadUnit_NO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadUnit_NO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadUnit_NO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadUnit_NO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadUnit_NO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadUnit_NO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadUnit_NO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadUnit_NO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadUnit_NO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadUnit_NO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadUnit_NO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadUnit_NO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadUnit_NO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadUnit_NO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadUnit_NO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadFile_TO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadFile_TO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadFile_TO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadFile_TO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadFile_TO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadFile_TO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadFile_TO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadFile_TO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadFile_TO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadFile_TO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadFile_TO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadFile_TO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadFile_TO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadFile_TO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadFile_TO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadFile_TO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadFile_TO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadFile_TO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadFile_TO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadFile_TO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadFile_TO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadFile_TO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadFile_TO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadFile_TO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadFile_TO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadUnit_TO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadUnit_TO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadUnit_TO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadUnit_TO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadUnit_TO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadUnit_TO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadUnit_TO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadUnit_TO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadUnit_TO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadUnit_TO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadUnit_TO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadUnit_TO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadUnit_TO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadUnit_TO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadUnit_TO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadUnit_TO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadUnit_TO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadUnit_TO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadUnit_TO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadUnit_TO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadUnit_TO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadUnit_TO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadUnit_TO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadUnit_TO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadUnit_TO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadFile_TO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadFile_TO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadFile_TO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadFile_TO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadFile_TO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadFile_TO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadFile_TO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadFile_TO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadFile_TO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadFile_TO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadFile_TO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadFile_TO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadFile_TO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadFile_TO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadFile_TO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadFile_TO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadFile_TO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadFile_TO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadFile_TO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadFile_TO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadFile_TO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadFile_TO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadFile_TO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadFile_TO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadFile_TO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableReadUnit_TO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableReadUnit_TO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableReadUnit_TO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableReadUnit_TO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableReadUnit_TO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableReadUnit_TO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableReadUnit_TO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableReadUnit_TO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableReadUnit_TO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableReadUnit_TO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableReadUnit_TO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableReadUnit_TO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableReadUnit_TO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableReadUnit_TO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableReadUnit_TO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableReadUnit_TO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableReadUnit_TO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableReadUnit_TO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableReadUnit_TO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableReadUnit_TO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableReadUnit_TO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableReadUnit_TO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableReadUnit_TO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableReadUnit_TO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableReadUnit_TO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getErrTableRead_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getErrTableWrite_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteFile_NO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteFile_NO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteFile_NO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteFile_NO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteFile_NO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteFile_NO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteFile_NO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteFile_NO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteFile_NO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteFile_NO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteFile_NO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteFile_NO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteFile_NO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteFile_NO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteFile_NO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteFile_NO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteFile_NO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteFile_NO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteFile_NO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteFile_NO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteFile_NO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteFile_NO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteFile_NO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteFile_NO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteFile_NO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteFile_NO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteFile_NO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteFile_NO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteFile_NO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteFile_NO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteFile_NO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteFile_NO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteFile_NO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteFile_NO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteFile_NO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteFile_NO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteFile_NO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteFile_NO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteFile_NO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteFile_NO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteFile_NO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteFile_NO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteFile_NO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteFile_NO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteFile_NO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteFile_NO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteFile_NO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteFile_NO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteFile_NO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteFile_NO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteUnit_NO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteFile_TO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteFile_TO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteFile_TO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteFile_TO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteFile_TO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteFile_TO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteFile_TO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteFile_TO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteFile_TO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteFile_TO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteFile_TO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteFile_TO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteFile_TO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteFile_TO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteFile_TO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteFile_TO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteFile_TO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteFile_TO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteFile_TO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteFile_TO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteFile_TO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteFile_TO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteFile_TO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteFile_TO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteFile_TO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteFile_TO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteFile_TO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteFile_TO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteFile_TO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteFile_TO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteFile_TO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteFile_TO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteFile_TO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteFile_TO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteFile_TO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteFile_TO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteFile_TO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteFile_TO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteFile_TO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteFile_TO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteFile_TO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteFile_TO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteFile_TO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteFile_TO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteFile_TO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteFile_TO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteFile_TO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteFile_TO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteFile_TO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteFile_TO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef File_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_IK5
        use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_IK4
        use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_IK3
        use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_IK2
        use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_IK1
        use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_LK5
        use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_LK4
        use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_LK3
        use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_LK2
        use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_LK1
        use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getErrTableWriteUnit_TO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getErrTableWrite_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFieldSep_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ID0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FDEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_ID0_FDEF_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_ID0_FDEF_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FDEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FCSV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_ID0_FCSV_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_ID0_FCSV_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FCSV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FFLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_ID0_FFLD_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_ID0_FFLD_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FFLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FDEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_ID0_FDEF_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_ID0_FDEF_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FDEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FCSV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_ID0_FCSV_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_ID0_FCSV_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FCSV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FFLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_ID0_FFLD_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_ID0_FFLD_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FFLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ID0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CD1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FDEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_CD1_FDEF_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_CD1_FDEF_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FDEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FCSV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_CD1_FCSV_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_CD1_FCSV_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FCSV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FFLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_CD1_FFLD_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_CD1_FFLD_XX_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FFLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FDEF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_CD1_FDEF_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_CD1_FDEF_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FDEF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FCSV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_CD1_FCSV_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_CD1_FCSV_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FCSV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FFLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define File_ENABLED 1
    module procedure getFieldSepFile_CD1_FFLD_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure getFieldSepUnit_CD1_FFLD_NF_SK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FFLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CD1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFieldSep_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define display_typer_ENABLED 1

#define File_ENABLED 1
    module procedure display_typer_file
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef File_ENABLED

#define Unit_ENABLED 1
    module procedure display_typer_unit
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure
#undef Unit_ENABLED

#undef display_typer_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !module procedure finalizeDisplay
    !    !>  \warning
    !    !>  Note that the finalize routine is called upon assignment and destroying of local copies of the derived type.<br>
    !    !>  The intel Fortran compiler 2021.5 **correctly** calls the final subroutine where needed, causing the opened
    !    !>  display unit to get closed. As such, the file remains closed after object construction which leads to runtime errors.
    !    !>  GNU Fortran compiler incorrectly does not call the `final` subroutine which leaves the file opened as desired,
    !    !>  but that is an incorrect non-standard-confoming behavior which should not be relied upon. The remedy as suggested by Jim in
    !    !>  [reported to the Intel team](https://community.intel.com/t5/Intel-Fortran-Compiler/COMPILER-BUG-duplicate-call-to-final-subroutine-for-object/m-p/1373155#M160832)
    !    !>  is to define an internal private component that keeps track of the open status of the file to write.
    !    if (.not. isPreconnected(self%unit)) then
    !        !inquire(unit = self%unit, opened = self%opened)
    !        if (self%opened) close(self%unit)
    !    end if
    !end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure skip_IK
        integer(IK) :: def_unit
        integer(IK) :: def_count
        if (present(sticky)) self%sticky = sticky
        if (present(unit )) then; def_unit  = unit  ; if (self%sticky) self%unit    = unit  ; else; def_unit    = self%unit     ; end if
        if (present(count)) then; def_count = count ; if (self%sticky) self%count   = count ; else; def_count   = self%count    ; end if
        if (0_IK < def_count) then
#if         MEXPRINT_ENABLED
            if (def_unit == output_unit) then
                call mexPrintf(repeat(NLC, def_count))
                return
            end if
#endif
            write(def_unit, "("//repeat("/", def_count)//")", advance = "no")
        end if
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getLenFieldMin_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getLenFieldMin_D0_IK5
!        use pm_kind, only: IKG => IK5
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getLenFieldMin_D0_IK4
!        use pm_kind, only: IKG => IK4
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getLenFieldMin_D0_IK3
!        use pm_kind, only: IKG => IK3
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getLenFieldMin_D0_IK2
!        use pm_kind, only: IKG => IK2
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getLenFieldMin_D0_IK1
!        use pm_kind, only: IKG => IK1
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define CK_ENABLED 1
!
!#if CK5_ENABLED
!    module procedure getLenFieldMin_D0_CK5
!        use pm_kind, only: CKG => CK5
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK4_ENABLED
!    module procedure getLenFieldMin_D0_CK4
!        use pm_kind, only: CKG => CK4
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK3_ENABLED
!    module procedure getLenFieldMin_D0_CK3
!        use pm_kind, only: CKG => CK3
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK2_ENABLED
!    module procedure getLenFieldMin_D0_CK2
!        use pm_kind, only: CKG => CK2
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if CK1_ENABLED
!    module procedure getLenFieldMin_D0_CK1
!        use pm_kind, only: CKG => CK1
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#undef CK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getLenFieldMin_D0_RK5
!        use pm_kind, only: RKG => RK5
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getLenFieldMin_D0_RK4
!        use pm_kind, only: RKG => RK4
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getLenFieldMin_D0_RK3
!        use pm_kind, only: RKG => RK3
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getLenFieldMin_D0_RK2
!        use pm_kind, only: RKG => RK2
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getLenFieldMin_D0_RK1
!        use pm_kind, only: RKG => RK1
!#include "pm_io@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getLenFieldMin_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED

#define getFormat_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    module procedure getFormat_D0_Def
#include "pm_io@routines.inc.F90"
    end procedure

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getFormat_D1_SK5
    use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getFormat_D1_SK4
    use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getFormat_D1_SK3
    use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getFormat_D1_SK2
    use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getFormat_D1_SK1
    use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getFormat_D1_IK5
    use pm_kind, only: IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getFormat_D1_IK4
    use pm_kind, only: IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getFormat_D1_IK3
    use pm_kind, only: IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getFormat_D1_IK2
    use pm_kind, only: IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getFormat_D1_IK1
    use pm_kind, only: IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getFormat_D1_LK5
    use pm_kind, only: LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getFormat_D1_LK4
    use pm_kind, only: LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getFormat_D1_LK3
    use pm_kind, only: LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getFormat_D1_LK2
    use pm_kind, only: LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getFormat_D1_LK1
    use pm_kind, only: LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getFormat_D1_CK5
    use pm_kind, only: CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getFormat_D1_CK4
    use pm_kind, only: CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getFormat_D1_CK3
    use pm_kind, only: CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getFormat_D1_CK2
    use pm_kind, only: CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getFormat_D1_CK1
    use pm_kind, only: CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getFormat_D1_RK5
    use pm_kind, only: RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getFormat_D1_RK4
    use pm_kind, only: RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getFormat_D1_RK3
    use pm_kind, only: RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getFormat_D1_RK2
    use pm_kind, only: RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getFormat_D1_RK1
    use pm_kind, only: RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFormat_ENABLED

#endif
!FORTRAN_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define wrap_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure wrap_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure wrap_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure wrap_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure wrap_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure wrap_D0_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef wrap_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define show_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Intrinsic Scalar
#define CN_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure show_D0_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D0_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D0_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D0_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D0_SK1
        use pm_kind, only: SKG => SK1 ! LCOV_EXCL_LINE
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure show_D0_IK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D0_IK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D0_IK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D0_IK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D0_IK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure show_D0_LK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D0_LK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D0_LK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D0_LK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D0_LK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure show_D0_CK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D0_CK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D0_CK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D0_CK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D0_CK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure show_D0_RK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D0_RK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D0_RK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D0_RK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D0_RK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure show_D1_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D1_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D1_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D1_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D1_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure show_D1_IK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D1_IK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D1_IK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D1_IK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D1_IK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure show_D1_LK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D1_LK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D1_LK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D1_LK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D1_LK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure show_D1_CK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D1_CK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D1_CK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D1_CK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D1_CK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure show_D1_RK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D1_RK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D1_RK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D1_RK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D1_RK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure show_D2_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D2_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D2_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D2_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D2_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure show_D2_IK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D2_IK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D2_IK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D2_IK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D2_IK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure show_D2_LK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D2_LK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D2_LK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D2_LK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D2_LK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure show_D2_CK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D2_CK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D2_CK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D2_CK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D2_CK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure show_D2_RK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D2_RK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D2_RK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D2_RK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D2_RK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure show_D3_SK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D3_SK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D3_SK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D3_SK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D3_SK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure show_D3_IK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D3_IK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D3_IK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D3_IK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D3_IK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure show_D3_LK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D3_LK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D3_LK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D3_LK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D3_LK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure show_D3_CK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D3_CK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D3_CK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D3_CK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D3_CK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure show_D3_RK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D3_RK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D3_RK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D3_RK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D3_RK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CN_ENABLED

    ! LCOV_EXCL_START
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D0_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D0_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D0_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D0_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D0_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D0_PSIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D0_PSIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D0_PSIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D0_PSIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D0_PSIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D0_PSLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D0_PSLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D0_PSLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D0_PSLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D0_PSLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D0_PSCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D0_PSCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D0_PSCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D0_PSCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D0_PSCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D0_PSRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D0_PSRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D0_PSRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D0_PSRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D0_PSRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D1_PSIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D1_PSIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D1_PSIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D1_PSIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D1_PSIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D1_PSLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D1_PSLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D1_PSLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D1_PSLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D1_PSLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D1_PSCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D1_PSCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D1_PSCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D1_PSCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D1_PSCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D1_PSRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D1_PSRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D1_PSRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D1_PSRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D1_PSRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D2_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D2_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D2_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D2_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D2_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D2_PSIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D2_PSIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D2_PSIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D2_PSIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D2_PSIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D2_PSLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D2_PSLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D2_PSLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D2_PSLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D2_PSLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D2_PSCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D2_PSCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D2_PSCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D2_PSCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D2_PSCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D2_PSRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D2_PSRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D2_PSRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D2_PSRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D2_PSRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D3_PSSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D3_PSSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D3_PSSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D3_PSSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D3_PSSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D3_PSIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D3_PSIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D3_PSIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D3_PSIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D3_PSIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D3_PSLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D3_PSLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D3_PSLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D3_PSLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D3_PSLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D3_PSCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D3_PSCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D3_PSCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D3_PSCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D3_PSCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D3_PSRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D3_PSRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D3_PSRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D3_PSRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D3_PSRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D0_PVSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D0_PVSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D0_PVSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D0_PVSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D0_PVSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D0_PVIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D0_PVIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D0_PVIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D0_PVIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D0_PVIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D0_PVLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D0_PVLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D0_PVLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D0_PVLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D0_PVLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D0_PVCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D0_PVCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D0_PVCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D0_PVCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D0_PVCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D0_PVRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D0_PVRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D0_PVRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D0_PVRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D0_PVRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D1_PVSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D1_PVSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D1_PVSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D1_PVSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D1_PVSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D1_PVIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D1_PVIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D1_PVIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D1_PVIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D1_PVIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D1_PVLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D1_PVLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D1_PVLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D1_PVLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D1_PVLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D1_PVCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D1_PVCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D1_PVCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D1_PVCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D1_PVCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D1_PVRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D1_PVRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D1_PVRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D1_PVRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D1_PVRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D2_PVSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D2_PVSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D2_PVSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D2_PVSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D2_PVSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D2_PVIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D2_PVIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D2_PVIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D2_PVIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D2_PVIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D2_PVLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D2_PVLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D2_PVLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D2_PVLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D2_PVLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D2_PVCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D2_PVCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D2_PVCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D2_PVCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D2_PVCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D2_PVRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D2_PVRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D2_PVRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D2_PVRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D2_PVRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D0_PMSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D0_PMSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D0_PMSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D0_PMSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D0_PMSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D0_PMIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D0_PMIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D0_PMIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D0_PMIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D0_PMIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D0_PMLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D0_PMLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D0_PMLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D0_PMLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D0_PMLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D0_PMCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D0_PMCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D0_PMCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D0_PMCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D0_PMCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D0_PMRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D0_PMRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D0_PMRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D0_PMRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D0_PMRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D1_PMSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D1_PMSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D1_PMSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D1_PMSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D1_PMSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D1_PMIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D1_PMIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D1_PMIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D1_PMIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D1_PMIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D1_PMLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D1_PMLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D1_PMLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D1_PMLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D1_PMLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D1_PMCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D1_PMCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D1_PMCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D1_PMCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D1_PMCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D1_PMRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D1_PMRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D1_PMRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D1_PMRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D1_PMRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if PDT_ENABLED
#if SK5_ENABLED
    module procedure show_D0_PCSK5
        use pm_kind, only: SKG => SK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure show_D0_PCSK4
        use pm_kind, only: SKG => SK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure show_D0_PCSK3
        use pm_kind, only: SKG => SK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure show_D0_PCSK2
        use pm_kind, only: SKG => SK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure show_D0_PCSK1
        use pm_kind, only: SKG => SK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if PDT_ENABLED
#if IK5_ENABLED
    module procedure show_D0_PCIK5
        use pm_kind, only: SKG => SK, IKG => IK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure show_D0_PCIK4
        use pm_kind, only: SKG => SK, IKG => IK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure show_D0_PCIK3
        use pm_kind, only: SKG => SK, IKG => IK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure show_D0_PCIK2
        use pm_kind, only: SKG => SK, IKG => IK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure show_D0_PCIK1
        use pm_kind, only: SKG => SK, IKG => IK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if PDT_ENABLED
#if LK5_ENABLED
    module procedure show_D0_PCLK5
        use pm_kind, only: SKG => SK, LKG => LK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure show_D0_PCLK4
        use pm_kind, only: SKG => SK, LKG => LK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure show_D0_PCLK3
        use pm_kind, only: SKG => SK, LKG => LK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure show_D0_PCLK2
        use pm_kind, only: SKG => SK, LKG => LK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure show_D0_PCLK1
        use pm_kind, only: SKG => SK, LKG => LK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if PDT_ENABLED
#if CK5_ENABLED
    module procedure show_D0_PCCK5
        use pm_kind, only: SKG => SK, CKG => CK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure show_D0_PCCK4
        use pm_kind, only: SKG => SK, CKG => CK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure show_D0_PCCK3
        use pm_kind, only: SKG => SK, CKG => CK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure show_D0_PCCK2
        use pm_kind, only: SKG => SK, CKG => CK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure show_D0_PCCK1
        use pm_kind, only: SKG => SK, CKG => CK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if PDT_ENABLED
#if RK5_ENABLED
    module procedure show_D0_PCRK5
        use pm_kind, only: SKG => SK, RKG => RK5
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure show_D0_PCRK4
        use pm_kind, only: SKG => SK, RKG => RK4
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure show_D0_PCRK3
        use pm_kind, only: SKG => SK, RKG => RK3
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure show_D0_PCRK2
        use pm_kind, only: SKG => SK, RKG => RK2
#include "pm_io@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure show_D0_PCRK1
        use pm_kind, only: SKG => SK, RKG => RK1
#include "pm_io@routines.inc.F90"
    end procedure
#endif
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PC_ENABLED

! LCOV_EXCL_STOP

! LCOV_EXCL_START

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Box of Scalar.
#define BS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D0_BSSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D0_BSIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D0_BSLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D0_BSCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D0_BSRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D1_BSSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D1_BSIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D1_BSLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D1_BSCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D1_BSRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D2_BSSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D2_BSIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D2_BSLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D2_BSCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D2_BSRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D3_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D3_BSSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D3_BSIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D3_BSLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D3_BSCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D3_BSRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D3_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D0_BVSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D0_BVIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D0_BVLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D0_BVCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D0_BVRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D1_BVSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D1_BVIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D1_BVLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D1_BVCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D1_BVRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D2_BVSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D2_BVIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D2_BVLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D2_BVCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D2_BVRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if FORTRAN_ENABLED

#define BM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D0_BMSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D0_BMIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D0_BMLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D0_BMCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D0_BMRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D1_BMSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D1_BMIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D1_BMLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D1_BMCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D1_BMRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    module procedure show_D0_BCSK
        use pm_kind, only: SKG => SK
#include "pm_io@routines.inc.F90"
    end procedure

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    module procedure show_D0_BCIK
        use pm_kind, only: SKG => SK, IKG => IK
#include "pm_io@routines.inc.F90"
    end procedure

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    module procedure show_D0_BCLK
        use pm_kind, only: SKG => SK, LKG => LK
#include "pm_io@routines.inc.F90"
    end procedure

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    module procedure show_D0_BCCK
        use pm_kind, only: SKG => SK, CKG => CK
#include "pm_io@routines.inc.F90"
    end procedure

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    module procedure show_D0_BCRK
        use pm_kind, only: SKG => SK, RKG => RK
#include "pm_io@routines.inc.F90"
    end procedure

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef BC_ENABLED

#endif
!FORTRAN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef show_ENABLED

! LCOV_EXCL_STOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! LCOV_EXCL_START

#if FORTRAN_ENABLED

#define dump_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure dump_D0
        use pm_kind ! must include all kinds
#include "pm_io@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure dump_D1
        use pm_kind ! must include all kinds
#include "pm_io@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure dump_D2
        use pm_kind ! must include all kinds
#include "pm_io@routines.inc.F90"
    end procedure

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef dump_ENABLED

#endif
!FORTRAN_ENABLED

! LCOV_EXCL_STOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines